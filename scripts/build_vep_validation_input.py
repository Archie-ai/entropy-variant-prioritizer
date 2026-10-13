#!/usr/bin/env python3
"""Create a deterministic 40-variant VEP validation VCF and manifest."""

import csv
import gzip
import hashlib
from collections import Counter
from pathlib import Path


INPUT = Path("data/processed/clinvar_benchmark_20260822.tsv.gz")
OUTPUT_DIR = Path("data/interim/vep_validation")
VCF_OUTPUT = OUTPUT_DIR / "clinvar_vep_validation_40.vcf"
MANIFEST_OUTPUT = OUTPUT_DIR / "clinvar_vep_validation_40_manifest.tsv"

EXPECTED_INPUT_SHA256 = (
    "2a74db014bdf5830b1af7c4cd296c9ea1c8ea8d9e91859d1302987b265e7176f"
)
EXPECTED_RECORDS = 65_222
PER_CLASS = 20
EXPECTED_LABELS = {"0": "benign", "1": "pathogenic"}
SPECIAL_CHROMOSOMES = ("X", "MT")

REQUIRED_COLUMNS = (
    "assembly",
    "clinvar_release",
    "chrom",
    "pos",
    "ref",
    "alt",
    "genomic_key",
    "clinvar_variation_id",
    "label",
    "label_name",
)


def sha256_file(path):
    """Calculate a file's SHA-256 checksum without loading it all at once."""
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def stable_rank(row):
    """Create a repeatable, row-order-independent selection value."""
    text = (
        f"{row['label']}|{row['genomic_key']}|"
        f"{row['clinvar_variation_id']}"
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def chromosome_sort_key(chrom):
    """Order chromosomes as 1-22, X, Y, MT, then unexpected values."""
    value = chrom.upper().removeprefix("CHR")

    if value.isdigit() and 1 <= int(value) <= 22:
        return (0, int(value))

    special_order = {"X": 23, "Y": 24, "MT": 25, "M": 25}

    if value in special_order:
        return (0, special_order[value])

    return (1, value)


def read_and_validate_benchmark():
    """Read the frozen benchmark and repeat its critical integrity checks."""
    if not INPUT.is_file():
        raise FileNotFoundError(f"Benchmark not found: {INPUT}")

    input_sha256 = sha256_file(INPUT)

    if input_sha256 != EXPECTED_INPUT_SHA256:
        raise ValueError(
            "Benchmark checksum mismatch:\n"
            f"  expected: {EXPECTED_INPUT_SHA256}\n"
            f"  observed: {input_sha256}"
        )

    with gzip.open(
        INPUT,
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("Benchmark header is missing")

        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    missing_columns = [
        column
        for column in REQUIRED_COLUMNS
        if column not in fieldnames
    ]

    if missing_columns:
        raise ValueError(f"Missing columns: {missing_columns}")

    if len(rows) != EXPECTED_RECORDS:
        raise ValueError(
            f"Expected {EXPECTED_RECORDS} records, "
            f"observed {len(rows)}"
        )

    if any(
        not row[column].strip()
        for row in rows
        for column in REQUIRED_COLUMNS
    ):
        raise ValueError("At least one required field is empty")

    if {row["assembly"] for row in rows} != {"GRCh38"}:
        raise ValueError("All records must use GRCh38")

    if any(
        not row["pos"].isdigit() or int(row["pos"]) < 1
        for row in rows
    ):
        raise ValueError("At least one genomic position is invalid")

    observed_labels = {
        (row["label"], row["label_name"])
        for row in rows
    }

    if observed_labels != set(EXPECTED_LABELS.items()):
        raise ValueError(
            f"Unexpected label mapping: {observed_labels}"
        )

    variation_ids = [
        row["clinvar_variation_id"]
        for row in rows
    ]
    genomic_keys = [
        row["genomic_key"]
        for row in rows
    ]

    if len(variation_ids) != len(set(variation_ids)):
        raise ValueError("ClinVar Variation IDs are not unique")

    if len(genomic_keys) != len(set(genomic_keys)):
        raise ValueError("Genomic keys are not unique")

    allowed_bases = {"A", "C", "G", "T"}

    for row in rows:
        ref = row["ref"].upper()
        alt = row["alt"].upper()

        if (
            len(ref) != 1
            or len(alt) != 1
            or ref not in allowed_bases
            or alt not in allowed_bases
            or ref == alt
        ):
            raise ValueError(
                f"Noncanonical SNV detected: "
                f"{row['genomic_key']}"
            )

    return rows, fieldnames, input_sha256


def select_validation_rows(rows):
    """Select 20 records per class using special cases and stable hashes."""
    selected = []

    for label in sorted(EXPECTED_LABELS):
        candidates = [
            row
            for row in rows
            if row["label"] == label
        ]

        ranked = sorted(
            candidates,
            key=lambda row: (
                stable_rank(row),
                row["genomic_key"],
            ),
        )

        class_selected = []
        chosen_keys = set()

        for chromosome in SPECIAL_CHROMOSOMES:
            special = [
                row
                for row in ranked
                if row["chrom"].upper() == chromosome
            ]

            if special:
                chosen = special[0].copy()
                chosen["selection_reason"] = (
                    f"special_chromosome_{chromosome}"
                )
                class_selected.append(chosen)
                chosen_keys.add(chosen["genomic_key"])

        for row in ranked:
            if len(class_selected) == PER_CLASS:
                break

            if row["genomic_key"] in chosen_keys:
                continue

            chosen = row.copy()
            chosen["selection_reason"] = "stable_hash_rank"
            class_selected.append(chosen)
            chosen_keys.add(chosen["genomic_key"])

        if len(class_selected) != PER_CLASS:
            raise ValueError(
                f"Could not select {PER_CLASS} records "
                f"for label {label}"
            )

        selected.extend(class_selected)

    for row in selected:
        row["selection_sha256"] = stable_rank(row)

    return sorted(
        selected,
        key=lambda row: (
            chromosome_sort_key(row["chrom"]),
            int(row["pos"]),
            row["ref"],
            row["alt"],
            row["clinvar_variation_id"],
        ),
    )


def write_outputs(rows, original_fieldnames):
    """Write a label-free VCF and a separate linkage manifest."""
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with VCF_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="\n",
    ) as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=entropy-variant-prioritizer\n")
        handle.write("##reference=GRCh38\n")
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )

        for row in rows:
            vep_input_id = (
                f"CV{row['clinvar_variation_id']}"
            )

            fields = (
                row["chrom"],
                row["pos"],
                vep_input_id,
                row["ref"].upper(),
                row["alt"].upper(),
                ".",
                "PASS",
                ".",
            )

            handle.write("\t".join(fields) + "\n")

    manifest_fields = [
        "vep_input_id",
        "vcf_order",
        "selection_reason",
        "selection_sha256",
        *original_fieldnames,
    ]

    with MANIFEST_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=manifest_fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()

        for vcf_order, row in enumerate(rows, start=1):
            output_row = dict(row)
            output_row["vep_input_id"] = (
                f"CV{row['clinvar_variation_id']}"
            )
            output_row["vcf_order"] = vcf_order
            writer.writerow(output_row)


def main():
    rows, fieldnames, input_sha256 = (
        read_and_validate_benchmark()
    )

    selected = select_validation_rows(rows)
    write_outputs(selected, fieldnames)

    print("VEP VALIDATION INPUT BUILD")
    print("--------------------------")
    print(f"Input records verified: {len(rows)}")
    print(f"Input SHA-256: {input_sha256}")
    print(
        f"Validation records written: {len(selected)}"
    )
    print(
        "Label counts: "
        f"{dict(Counter(row['label_name'] for row in selected))}"
    )
    print(
        "Chromosome counts: "
        f"{dict(Counter(row['chrom'] for row in selected))}"
    )
    print(
        "Selection reasons: "
        f"{dict(Counter(row['selection_reason'] for row in selected))}"
    )
    print(f"VCF: {VCF_OUTPUT}")
    print(
        f"VCF SHA-256: {sha256_file(VCF_OUTPUT)}"
    )
    print(f"Manifest: {MANIFEST_OUTPUT}")
    print(
        "Manifest SHA-256: "
        f"{sha256_file(MANIFEST_OUTPUT)}"
    )


if __name__ == "__main__":
    main()

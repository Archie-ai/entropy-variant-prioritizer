import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


BENCHMARK_FILE = Path(
    "data/processed/clinvar_benchmark_20260822.tsv.gz"
)
SCHEMA_FILE = Path(
    "data/processed/clinvar_benchmark_20260822.schema.json"
)
CHECKSUM_FILE = Path(
    "data/processed/clinvar_benchmark_20260822.sha256"
)
VALIDATION_REPORT = Path(
    "data/processed/"
    "clinvar_benchmark_20260822.validation.json"
)

ASSEMBLY = "GRCh38"
CLINVAR_RELEASE = "2026-08-22"

EXPECTED_SHA256 = (
    "2a74db014bdf5830b1af7c4cd296c9ea1"
    "c8ea8d9e91859d1302987b265e7176f"
)

CANONICAL_BASES = {"A", "C", "G", "T"}

PATHOGENIC_CLASSES = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
}

BENIGN_CLASSES = {
    "Benign",
    "Likely_benign",
    "Benign/Likely_benign",
}

REVIEW_STARS = {
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
}

EXPECTED_COLUMNS = [
    "assembly",
    "clinvar_release",
    "chrom",
    "pos",
    "ref",
    "alt",
    "genomic_key",
    "clinvar_variation_id",
    "clinvar_allele_id",
    "dbsnp_rs_id",
    "gene_info",
    "gene_symbols",
    "gene_ids",
    "molecular_consequence",
    "clinical_hgvs",
    "clinical_significance",
    "label",
    "label_name",
    "review_status",
    "review_stars",
    "variant_type",
]

EXPECTED_CLASS_COUNTS = {
    "pathogenic": 24904,
    "benign": 40318,
}

EXPECTED_REVIEW_COUNTS = {
    (2, "pathogenic"): 21535,
    (2, "benign"): 38813,
    (3, "pathogenic"): 3361,
    (3, "benign"): 1505,
    (4, "pathogenic"): 8,
    (4, "benign"): 0,
}

EXPECTED_RECORD_COUNT = 65222

OPTIONAL_ANNOTATIONS = [
    "clinvar_allele_id",
    "dbsnp_rs_id",
    "gene_info",
    "clinical_hgvs",
]

MAXIMUM_ERROR_EXAMPLES = 25


def calculate_sha256(file_path):
    """Calculate the SHA-256 checksum of a file."""

    digest = hashlib.sha256()

    with file_path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def extract_gene_fields(gene_info):
    """Independently reconstruct gene symbols and IDs from GENEINFO."""

    if not gene_info:
        return "", ""

    gene_symbols = []
    gene_ids = []

    for gene_entry in gene_info.split("|"):
        gene_entry = gene_entry.strip()

        if not gene_entry:
            continue

        if ":" in gene_entry:
            gene_symbol, gene_id = gene_entry.rsplit(":", 1)
        else:
            gene_symbol = gene_entry
            gene_id = ""

        gene_symbols.append(gene_symbol)
        gene_ids.append(gene_id)

    return "|".join(gene_symbols), "|".join(gene_ids)


def expected_classification(clinical_significance):
    """Determine the expected label from the original CLNSIG value."""

    if clinical_significance in PATHOGENIC_CLASSES:
        return 1, "pathogenic"

    if clinical_significance in BENIGN_CLASSES:
        return 0, "benign"

    return None


def write_json_atomically(file_path, data):
    """Write JSON through a temporary file before replacing the target."""

    temporary_path = file_path.with_name(file_path.name + ".tmp")

    with temporary_path.open(
        "w",
        encoding="utf-8",
        newline="\n",
    ) as file:
        json.dump(
            data,
            file,
            indent=2,
            ensure_ascii=False,
        )
        file.write("\n")

    temporary_path.replace(file_path)


def main():
    required_files = [
        BENCHMARK_FILE,
        SCHEMA_FILE,
        CHECKSUM_FILE,
    ]

    for required_file in required_files:
        if not required_file.exists():
            raise FileNotFoundError(
                f"Required file was not found: {required_file}"
            )

    error_count = 0
    error_examples = []

    def add_error(message):
        nonlocal error_count

        error_count += 1

        if len(error_examples) < MAXIMUM_ERROR_EXAMPLES:
            error_examples.append(message)

    actual_sha256 = calculate_sha256(BENCHMARK_FILE)

    if actual_sha256 != EXPECTED_SHA256:
        add_error(
            "Benchmark SHA-256 does not match the frozen expected "
            f"value. Expected {EXPECTED_SHA256}; "
            f"observed {actual_sha256}."
        )

    checksum_text = CHECKSUM_FILE.read_text(
        encoding="utf-8"
    ).strip()

    checksum_parts = checksum_text.split()

    if len(checksum_parts) != 2:
        add_error(
            "Checksum manifest does not contain exactly a checksum "
            "and filename."
        )
        manifest_sha256 = ""
        manifest_filename = ""
    else:
        manifest_sha256 = checksum_parts[0]
        manifest_filename = checksum_parts[1]

        if manifest_sha256 != actual_sha256:
            add_error(
                "Checksum manifest SHA-256 does not match the "
                "benchmark file."
            )

        if manifest_filename != BENCHMARK_FILE.name:
            add_error(
                "Checksum manifest filename does not match the "
                "benchmark filename."
            )

    try:
        schema = json.loads(
            SCHEMA_FILE.read_text(encoding="utf-8")
        )
    except json.JSONDecodeError as error:
        add_error(f"Schema JSON could not be parsed: {error}")
        schema = {}

    if schema.get("output_file") != BENCHMARK_FILE.name:
        add_error(
            "Schema output filename does not match the benchmark "
            "filename."
        )

    if schema.get("output_sha256") != actual_sha256:
        add_error(
            "Schema SHA-256 does not match the benchmark file."
        )

    if schema.get("assembly") != ASSEMBLY:
        add_error(
            f"Schema assembly is not {ASSEMBLY}."
        )

    if schema.get("source_release") != CLINVAR_RELEASE:
        add_error(
            "Schema ClinVar release does not match the frozen "
            f"release {CLINVAR_RELEASE}."
        )

    schema_columns = []

    for column_definition in schema.get("columns", []):
        if isinstance(column_definition, dict):
            schema_columns.append(
                column_definition.get("name")
            )

    if schema_columns != EXPECTED_COLUMNS:
        add_error(
            "Schema column order does not match the expected "
            "benchmark column order."
        )

    total_records = 0
    class_counts = Counter()
    review_counts = Counter()
    missing_annotation_counts = Counter()

    genomic_keys = set()
    variation_ids = set()
    allele_ids = set()

    duplicate_genomic_keys = 0
    duplicate_variation_ids = 0
    duplicate_allele_id_excess = 0

    with gzip.open(
        BENCHMARK_FILE,
        "rt",
        encoding="utf-8",
        newline="",
    ) as file:
        reader = csv.DictReader(
            file,
            delimiter="\t",
        )

        observed_columns = reader.fieldnames or []

        if observed_columns != EXPECTED_COLUMNS:
            add_error(
                "TSV header does not exactly match the expected "
                "column names and order."
            )

        for source_line_number, row in enumerate(
            reader,
            start=2,
        ):
            total_records += 1

            if None in row:
                add_error(
                    f"Line {source_line_number}: extra unnamed "
                    "columns were detected."
                )

            missing_structural_values = [
                column_name
                for column_name in EXPECTED_COLUMNS
                if row.get(column_name) is None
            ]

            if missing_structural_values:
                add_error(
                    f"Line {source_line_number}: missing columns "
                    f"{missing_structural_values}."
                )
                continue

            assembly = row["assembly"]
            release = row["clinvar_release"]
            chrom = row["chrom"]
            pos_text = row["pos"]
            ref = row["ref"]
            alt = row["alt"]
            genomic_key = row["genomic_key"]
            variation_id = row["clinvar_variation_id"]
            allele_id = row["clinvar_allele_id"]
            gene_info = row["gene_info"]
            gene_symbols = row["gene_symbols"]
            gene_ids = row["gene_ids"]
            molecular_consequence = row[
                "molecular_consequence"
            ]
            clinical_significance = row[
                "clinical_significance"
            ]
            label_text = row["label"]
            label_name = row["label_name"]
            review_status = row["review_status"]
            review_stars_text = row["review_stars"]
            variant_type = row["variant_type"]

            if assembly != ASSEMBLY:
                add_error(
                    f"Line {source_line_number}: unexpected "
                    f"assembly {assembly!r}."
                )

            if release != CLINVAR_RELEASE:
                add_error(
                    f"Line {source_line_number}: unexpected "
                    f"ClinVar release {release!r}."
                )

            if not chrom:
                add_error(
                    f"Line {source_line_number}: chromosome is "
                    "missing."
                )

            try:
                pos = int(pos_text)

                if pos <= 0:
                    add_error(
                        f"Line {source_line_number}: genomic "
                        "position is not positive."
                    )

            except ValueError:
                pos = None
                add_error(
                    f"Line {source_line_number}: position "
                    f"{pos_text!r} is not an integer."
                )

            if ref not in CANONICAL_BASES:
                add_error(
                    f"Line {source_line_number}: noncanonical "
                    f"REF allele {ref!r}."
                )

            if alt not in CANONICAL_BASES:
                add_error(
                    f"Line {source_line_number}: noncanonical "
                    f"ALT allele {alt!r}."
                )

            if ref == alt:
                add_error(
                    f"Line {source_line_number}: REF and ALT "
                    "are identical."
                )

            reconstructed_key = (
                f"{chrom}:{pos_text}:{ref}:{alt}"
            )

            if genomic_key != reconstructed_key:
                add_error(
                    f"Line {source_line_number}: genomic key "
                    f"{genomic_key!r} does not match "
                    f"{reconstructed_key!r}."
                )

            if genomic_key in genomic_keys:
                duplicate_genomic_keys += 1
                add_error(
                    f"Line {source_line_number}: duplicate "
                    f"genomic key {genomic_key!r}."
                )
            else:
                genomic_keys.add(genomic_key)

            if not variation_id or variation_id == ".":
                add_error(
                    f"Line {source_line_number}: ClinVar "
                    "Variation ID is missing."
                )
            elif variation_id in variation_ids:
                duplicate_variation_ids += 1
                add_error(
                    f"Line {source_line_number}: duplicate "
                    f"ClinVar Variation ID {variation_id!r}."
                )
            else:
                variation_ids.add(variation_id)

            if allele_id:
                if allele_id in allele_ids:
                    duplicate_allele_id_excess += 1
                else:
                    allele_ids.add(allele_id)

            expected_label = expected_classification(
                clinical_significance
            )

            if expected_label is None:
                add_error(
                    f"Line {source_line_number}: unaccepted "
                    f"clinical significance "
                    f"{clinical_significance!r}."
                )
            else:
                expected_label_number, expected_label_name = (
                    expected_label
                )

                try:
                    observed_label_number = int(label_text)
                except ValueError:
                    observed_label_number = None
                    add_error(
                        f"Line {source_line_number}: label "
                        f"{label_text!r} is not an integer."
                    )

                if (
                    observed_label_number
                    != expected_label_number
                ):
                    add_error(
                        f"Line {source_line_number}: label does "
                        "not agree with clinical significance."
                    )

                if label_name != expected_label_name:
                    add_error(
                        f"Line {source_line_number}: label name "
                        "does not agree with clinical "
                        "significance."
                    )

            if label_name not in {"pathogenic", "benign"}:
                add_error(
                    f"Line {source_line_number}: invalid label "
                    f"name {label_name!r}."
                )
            else:
                class_counts[label_name] += 1

            expected_stars = REVIEW_STARS.get(review_status)

            if expected_stars is None:
                add_error(
                    f"Line {source_line_number}: unaccepted "
                    f"review status {review_status!r}."
                )

            try:
                observed_stars = int(review_stars_text)
            except ValueError:
                observed_stars = None
                add_error(
                    f"Line {source_line_number}: review-star "
                    f"value {review_stars_text!r} is not an "
                    "integer."
                )

            if observed_stars != expected_stars:
                add_error(
                    f"Line {source_line_number}: review stars "
                    "do not agree with review status."
                )

            if (
                observed_stars in {2, 3, 4}
                and label_name in {"pathogenic", "benign"}
            ):
                review_counts[
                    (observed_stars, label_name)
                ] += 1

            if variant_type != "single_nucleotide_variant":
                add_error(
                    f"Line {source_line_number}: unexpected "
                    f"variant type {variant_type!r}."
                )

            if (
                "missense_variant"
                not in molecular_consequence
            ):
                add_error(
                    f"Line {source_line_number}: molecular "
                    "consequence does not contain "
                    "missense_variant."
                )

            reconstructed_symbols, reconstructed_ids = (
                extract_gene_fields(gene_info)
            )

            if gene_symbols != reconstructed_symbols:
                add_error(
                    f"Line {source_line_number}: gene symbols "
                    "do not agree with GENEINFO."
                )

            if gene_ids != reconstructed_ids:
                add_error(
                    f"Line {source_line_number}: gene IDs do "
                    "not agree with GENEINFO."
                )

            for annotation_name in OPTIONAL_ANNOTATIONS:
                if not row[annotation_name]:
                    missing_annotation_counts[
                        annotation_name
                    ] += 1

    if total_records != EXPECTED_RECORD_COUNT:
        add_error(
            "Exported record count does not match the expected "
            f"{EXPECTED_RECORD_COUNT}. Observed {total_records}."
        )

    if len(genomic_keys) != EXPECTED_RECORD_COUNT:
        add_error(
            "Unique genomic-key count does not equal the expected "
            f"{EXPECTED_RECORD_COUNT}. Observed "
            f"{len(genomic_keys)}."
        )

    if len(variation_ids) != EXPECTED_RECORD_COUNT:
        add_error(
            "Unique ClinVar Variation-ID count does not equal "
            f"{EXPECTED_RECORD_COUNT}. Observed "
            f"{len(variation_ids)}."
        )

    for label_name, expected_count in (
        EXPECTED_CLASS_COUNTS.items()
    ):
        observed_count = class_counts[label_name]

        if observed_count != expected_count:
            add_error(
                f"{label_name} count: expected "
                f"{expected_count}, observed {observed_count}."
            )

    for review_key, expected_count in (
        EXPECTED_REVIEW_COUNTS.items()
    ):
        observed_count = review_counts[review_key]

        if observed_count != expected_count:
            stars, label_name = review_key

            add_error(
                f"{stars}-star {label_name} count: expected "
                f"{expected_count}, observed {observed_count}."
            )

    if schema.get("record_count") != total_records:
        add_error(
            "Schema record count does not match the TSV record "
            "count."
        )

    schema_class_counts = schema.get(
        "class_counts",
        {},
    )

    for label_name in ("pathogenic", "benign"):
        if (
            schema_class_counts.get(label_name)
            != class_counts[label_name]
        ):
            add_error(
                f"Schema {label_name} count does not match the "
                "TSV."
            )

    schema_review_counts = schema.get(
        "review_status_counts",
        {},
    )

    for stars in (2, 3, 4):
        schema_star_group = schema_review_counts.get(
            f"{stars}_star",
            {},
        )

        for label_name in ("pathogenic", "benign"):
            if (
                schema_star_group.get(label_name)
                != review_counts[(stars, label_name)]
            ):
                add_error(
                    f"Schema {stars}-star {label_name} count "
                    "does not match the TSV."
                )

    schema_missing_counts = schema.get(
        "missing_optional_annotations",
        {},
    )

    for annotation_name in OPTIONAL_ANNOTATIONS:
        schema_missing_count = schema_missing_counts.get(
            annotation_name,
            0,
        )
        observed_missing_count = missing_annotation_counts[
            annotation_name
        ]

        if schema_missing_count != observed_missing_count:
            add_error(
                f"Schema missing count for {annotation_name} "
                "does not match the TSV."
            )

    status = "PASS" if error_count == 0 else "FAIL"

    report = {
        "status": status,
        "validated_file": BENCHMARK_FILE.as_posix(),
        "compressed_size_bytes": BENCHMARK_FILE.stat().st_size,
        "sha256": actual_sha256,
        "record_count": total_records,
        "unique_genomic_keys": len(genomic_keys),
        "duplicate_genomic_keys": duplicate_genomic_keys,
        "unique_clinvar_variation_ids": len(variation_ids),
        "duplicate_clinvar_variation_ids": (
            duplicate_variation_ids
        ),
        "unique_nonmissing_clinvar_allele_ids": len(
            allele_ids
        ),
        "repeated_clinvar_allele_id_excess": (
            duplicate_allele_id_excess
        ),
        "class_counts": {
            "pathogenic": class_counts["pathogenic"],
            "benign": class_counts["benign"],
        },
        "review_status_counts": {
            "2_star": {
                "pathogenic": review_counts[
                    (2, "pathogenic")
                ],
                "benign": review_counts[(2, "benign")],
            },
            "3_star": {
                "pathogenic": review_counts[
                    (3, "pathogenic")
                ],
                "benign": review_counts[(3, "benign")],
            },
            "4_star": {
                "pathogenic": review_counts[
                    (4, "pathogenic")
                ],
                "benign": review_counts[(4, "benign")],
            },
        },
        "missing_optional_annotations": {
            annotation_name: missing_annotation_counts[
                annotation_name
            ]
            for annotation_name in OPTIONAL_ANNOTATIONS
        },
        "validation_error_count": error_count,
        "validation_error_examples": error_examples,
    }

    write_json_atomically(
        VALIDATION_REPORT,
        report,
    )

    print("INDEPENDENT BENCHMARK EXPORT VALIDATION")
    print("---------------------------------------")
    print("Validated file:", BENCHMARK_FILE)
    print("Validation report:", VALIDATION_REPORT)
    print("Compressed size in bytes:", BENCHMARK_FILE.stat().st_size)
    print("SHA-256:", actual_sha256)
    print()

    print("RECORD AND IDENTIFIER INTEGRITY")
    print("-------------------------------")
    print("Total records:", total_records)
    print("Unique genomic keys:", len(genomic_keys))
    print(
        "Duplicate genomic-key excess:",
        duplicate_genomic_keys,
    )
    print(
        "Unique ClinVar Variation IDs:",
        len(variation_ids),
    )
    print(
        "Duplicate ClinVar Variation-ID excess:",
        duplicate_variation_ids,
    )
    print(
        "Unique nonmissing ClinVar Allele IDs:",
        len(allele_ids),
    )
    print(
        "Repeated ClinVar Allele-ID excess:",
        duplicate_allele_id_excess,
    )
    print()

    print("CLASS COUNTS")
    print("------------")
    print(
        "Pathogenic / likely pathogenic:",
        class_counts["pathogenic"],
    )
    print(
        "Benign / likely benign:",
        class_counts["benign"],
    )
    print(
        "Class-count verification:",
        sum(class_counts.values()),
    )
    print()

    print("OPTIONAL ANNOTATION MISSING COUNTS")
    print("----------------------------------")

    for annotation_name in OPTIONAL_ANNOTATIONS:
        print(
            f"{annotation_name}:",
            missing_annotation_counts[annotation_name],
        )

    print()

    if error_count == 0:
        print("FINAL VALIDATION")
        print("----------------")
        print(
            "PASS: Checksum, schema, record structure, labels, "
            "review levels, identifiers, and uniqueness were "
            "independently verified for all 65,222 records."
        )
    else:
        print("FINAL VALIDATION")
        print("----------------")
        print(
            f"FAIL: {error_count} validation problem(s) were "
            "detected."
        )

        print()
        print("ERROR EXAMPLES")
        print("--------------")

        for error_message in error_examples:
            print("-", error_message)

        raise SystemExit(1)


if __name__ == "__main__":
    main()
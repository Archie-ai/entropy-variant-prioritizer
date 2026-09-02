import csv
import gzip
import hashlib
import io
import json
from collections import Counter
from pathlib import Path


INPUT_VCF = Path("data/raw/clinvar_20260822.vcf.gz")
OUTPUT_DIRECTORY = Path("data/processed")

OUTPUT_TSV = OUTPUT_DIRECTORY / "clinvar_benchmark_20260822.tsv.gz"
OUTPUT_SCHEMA = OUTPUT_DIRECTORY / "clinvar_benchmark_20260822.schema.json"
OUTPUT_CHECKSUM = OUTPUT_DIRECTORY / "clinvar_benchmark_20260822.sha256"

ASSEMBLY = "GRCh38"
CLINVAR_RELEASE = "2026-08-22"

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

CANONICAL_BASES = {"A", "C", "G", "T"}

EXPECTED_COUNTS = {
    "total_variants": 4467990,
    "missense_variants": 2536550,
    "missense_snvs": 2527633,
    "one_character_alleles": 2527633,
    "canonical_snvs": 2527597,
    "accepted_clinical_significance": 223891,
    "strict_records": 65222,
    "pathogenic_records": 24904,
    "benign_records": 40318,
}

EXPECTED_REVIEW_COUNTS = {
    (2, "pathogenic"): 21535,
    (2, "benign"): 38813,
    (3, "pathogenic"): 3361,
    (3, "benign"): 1505,
    (4, "pathogenic"): 8,
    (4, "benign"): 0,
}

OUTPUT_COLUMNS = [
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

COLUMN_SCHEMA = [
    {
        "name": "assembly",
        "type": "string",
        "description": "Reference genome assembly.",
    },
    {
        "name": "clinvar_release",
        "type": "string",
        "description": "Frozen ClinVar release date.",
    },
    {
        "name": "chrom",
        "type": "string",
        "description": "Chromosome from the ClinVar VCF.",
    },
    {
        "name": "pos",
        "type": "integer",
        "description": "One-based genomic position.",
    },
    {
        "name": "ref",
        "type": "string",
        "description": "Canonical reference nucleotide.",
    },
    {
        "name": "alt",
        "type": "string",
        "description": "Canonical alternate nucleotide.",
    },
    {
        "name": "genomic_key",
        "type": "string",
        "description": "Unique CHROM:POS:REF:ALT identifier.",
    },
    {
        "name": "clinvar_variation_id",
        "type": "string",
        "description": "ClinVar Variation ID from VCF column 3.",
    },
    {
        "name": "clinvar_allele_id",
        "type": "string",
        "description": "ClinVar Allele ID from INFO/ALLELEID.",
    },
    {
        "name": "dbsnp_rs_id",
        "type": "string",
        "description": "Raw dbSNP identifier value from INFO/RS; blank if absent.",
    },
    {
        "name": "gene_info",
        "type": "string",
        "description": "Unmodified INFO/GENEINFO value.",
    },
    {
        "name": "gene_symbols",
        "type": "string",
        "description": "Pipe-separated gene symbols extracted from GENEINFO.",
    },
    {
        "name": "gene_ids",
        "type": "string",
        "description": "Pipe-separated NCBI Gene IDs extracted from GENEINFO.",
    },
    {
        "name": "molecular_consequence",
        "type": "string",
        "description": "Unmodified INFO/MC molecular consequence annotation.",
    },
    {
        "name": "clinical_hgvs",
        "type": "string",
        "description": "Unmodified INFO/CLNHGVS value.",
    },
    {
        "name": "clinical_significance",
        "type": "string",
        "description": "Original accepted INFO/CLNSIG value.",
    },
    {
        "name": "label",
        "type": "integer",
        "description": "Binary classification label: 1 pathogenic, 0 benign.",
    },
    {
        "name": "label_name",
        "type": "string",
        "description": "Human-readable binary class name.",
    },
    {
        "name": "review_status",
        "type": "string",
        "description": "Unmodified INFO/CLNREVSTAT value.",
    },
    {
        "name": "review_stars",
        "type": "integer",
        "description": "ClinVar review confidence level: 2, 3, or 4 stars.",
    },
    {
        "name": "variant_type",
        "type": "string",
        "description": "ClinVar variant type from INFO/CLNVC.",
    },
]


def parse_info(info_text):
    """Parse a VCF INFO column into a dictionary."""

    info_dict = {}

    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = True

    return info_dict


def clean_optional_value(value):
    """Represent missing optional VCF values as an empty string."""

    if value is None or value is True or value in {"", "."}:
        return ""

    return str(value)


def extract_gene_fields(gene_info):
    """
    Extract gene symbols and NCBI Gene IDs from GENEINFO.

    Multiple annotations remain pipe-separated and in their original order.
    The original GENEINFO value is also retained separately in the output.
    """

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


def classify_clinical_significance(clinical_significance):
    """Return the binary label and label name for an accepted CLNSIG value."""

    if clinical_significance in PATHOGENIC_CLASSES:
        return 1, "pathogenic"

    if clinical_significance in BENIGN_CLASSES:
        return 0, "benign"

    return None


def calculate_sha256(file_path):
    """Calculate the SHA-256 checksum of a file."""

    digest = hashlib.sha256()

    with file_path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def write_text_atomically(file_path, text):
    """Write a text file through a temporary file before replacing the target."""

    temporary_path = file_path.with_name(file_path.name + ".tmp")

    with temporary_path.open(
        "w",
        encoding="utf-8",
        newline="\n",
    ) as file:
        file.write(text)

    temporary_path.replace(file_path)


def verify_counts(observed_counts, review_counts):
    """Stop the export if any frozen benchmark count has changed."""

    errors = []

    for count_name, expected_value in EXPECTED_COUNTS.items():
        observed_value = observed_counts[count_name]

        if observed_value != expected_value:
            errors.append(
                f"{count_name}: expected {expected_value}, "
                f"observed {observed_value}"
            )

    for review_key, expected_value in EXPECTED_REVIEW_COUNTS.items():
        observed_value = review_counts[review_key]

        if observed_value != expected_value:
            stars, label_name = review_key
            errors.append(
                f"{stars}-star {label_name}: expected {expected_value}, "
                f"observed {observed_value}"
            )

    if errors:
        error_message = "\n".join(f"- {error}" for error in errors)

        raise RuntimeError(
            "Benchmark export verification failed:\n" + error_message
        )


def build_schema(
    observed_counts,
    review_counts,
    annotation_missing_counts,
    sha256_checksum,
):
    """Construct the machine-readable benchmark schema and provenance record."""

    review_breakdown = {
        "2_star": {
            "pathogenic": review_counts[(2, "pathogenic")],
            "benign": review_counts[(2, "benign")],
        },
        "3_star": {
            "pathogenic": review_counts[(3, "pathogenic")],
            "benign": review_counts[(3, "benign")],
        },
        "4_star": {
            "pathogenic": review_counts[(4, "pathogenic")],
            "benign": review_counts[(4, "benign")],
        },
    }

    return {
        "dataset_name": "ClinVar strict missense SNV benchmark",
        "description": (
            "Canonical GRCh38 missense SNVs with pathogenic/likely "
            "pathogenic or benign/likely benign ClinVar classifications "
            "and 2-4-star review status."
        ),
        "source_vcf": INPUT_VCF.as_posix(),
        "source_release": CLINVAR_RELEASE,
        "assembly": ASSEMBLY,
        "output_file": OUTPUT_TSV.name,
        "output_sha256": sha256_checksum,
        "record_count": observed_counts["strict_records"],
        "class_counts": {
            "pathogenic": observed_counts["pathogenic_records"],
            "benign": observed_counts["benign_records"],
        },
        "review_status_counts": review_breakdown,
        "missing_optional_annotations": dict(annotation_missing_counts),
        "filter_definition": [
            "INFO/MC contains missense_variant",
            "INFO/CLNVC equals single_nucleotide_variant",
            "REF and ALT each contain exactly one character",
            "REF and ALT are both members of A/C/G/T",
            "INFO/CLNSIG belongs to an accepted binary class",
            "INFO/CLNREVSTAT represents a 2-star, 3-star, or 4-star review",
        ],
        "columns": COLUMN_SCHEMA,
    }


def main():
    if not INPUT_VCF.exists():
        raise FileNotFoundError(
            f"Input ClinVar VCF was not found: {INPUT_VCF}"
        )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    temporary_output = OUTPUT_TSV.with_name(OUTPUT_TSV.name + ".tmp")

    observed_counts = Counter()
    review_counts = Counter()
    annotation_missing_counts = Counter()
    genomic_keys = set()

    try:
        with gzip.open(
            INPUT_VCF,
            "rt",
            encoding="utf-8",
        ) as input_file:
            with temporary_output.open("wb") as raw_output:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    fileobj=raw_output,
                    mtime=0,
                ) as gzip_output:
                    with io.TextIOWrapper(
                        gzip_output,
                        encoding="utf-8",
                        newline="",
                    ) as text_output:
                        writer = csv.DictWriter(
                            text_output,
                            fieldnames=OUTPUT_COLUMNS,
                            delimiter="\t",
                            lineterminator="\n",
                            quoting=csv.QUOTE_MINIMAL,
                        )

                        writer.writeheader()

                        for line_number, line in enumerate(
                            input_file,
                            start=1,
                        ):
                            if line.startswith("#"):
                                continue

                            observed_counts["total_variants"] += 1

                            fields = line.rstrip("\n").split("\t")

                            if len(fields) < 8:
                                raise ValueError(
                                    "Malformed VCF record at source line "
                                    f"{line_number}: fewer than 8 columns"
                                )

                            chrom = fields[0]
                            pos = fields[1]
                            variation_id = fields[2]
                            ref = fields[3].upper()
                            alt = fields[4].upper()
                            info_dict = parse_info(fields[7])

                            molecular_consequence = str(
                                info_dict.get("MC", "")
                            )
                            variant_type = str(
                                info_dict.get("CLNVC", "")
                            )
                            clinical_significance = str(
                                info_dict.get("CLNSIG", "")
                            )
                            review_status = str(
                                info_dict.get("CLNREVSTAT", "")
                            )

                            if (
                                "missense_variant"
                                not in molecular_consequence
                            ):
                                continue

                            observed_counts["missense_variants"] += 1

                            if (
                                variant_type
                                != "single_nucleotide_variant"
                            ):
                                continue

                            observed_counts["missense_snvs"] += 1

                            if len(ref) != 1 or len(alt) != 1:
                                continue

                            observed_counts["one_character_alleles"] += 1

                            if (
                                ref not in CANONICAL_BASES
                                or alt not in CANONICAL_BASES
                            ):
                                continue

                            observed_counts["canonical_snvs"] += 1

                            classification = (
                                classify_clinical_significance(
                                    clinical_significance
                                )
                            )

                            if classification is None:
                                continue

                            observed_counts[
                                "accepted_clinical_significance"
                            ] += 1

                            if review_status not in REVIEW_STARS:
                                continue

                            label, label_name = classification
                            review_stars = REVIEW_STARS[review_status]

                            if variation_id in {"", "."}:
                                raise ValueError(
                                    "Strict record is missing a ClinVar "
                                    f"Variation ID at {chrom}:{pos}:{ref}:{alt}"
                                )

                            genomic_key = (
                                f"{chrom}:{pos}:{ref}:{alt}"
                            )

                            if genomic_key in genomic_keys:
                                raise ValueError(
                                    "Duplicate strict genomic key found: "
                                    f"{genomic_key}"
                                )

                            genomic_keys.add(genomic_key)

                            allele_id = clean_optional_value(
                                info_dict.get("ALLELEID")
                            )
                            rs_id = clean_optional_value(
                                info_dict.get("RS")
                            )
                            gene_info = clean_optional_value(
                                info_dict.get("GENEINFO")
                            )
                            clinical_hgvs = clean_optional_value(
                                info_dict.get("CLNHGVS")
                            )

                            gene_symbols, gene_ids = extract_gene_fields(
                                gene_info
                            )

                            optional_annotations = {
                                "clinvar_allele_id": allele_id,
                                "dbsnp_rs_id": rs_id,
                                "gene_info": gene_info,
                                "clinical_hgvs": clinical_hgvs,
                            }

                            for (
                                annotation_name,
                                annotation_value,
                            ) in optional_annotations.items():
                                if not annotation_value:
                                    annotation_missing_counts[
                                        annotation_name
                                    ] += 1

                            writer.writerow(
                                {
                                    "assembly": ASSEMBLY,
                                    "clinvar_release": CLINVAR_RELEASE,
                                    "chrom": chrom,
                                    "pos": pos,
                                    "ref": ref,
                                    "alt": alt,
                                    "genomic_key": genomic_key,
                                    "clinvar_variation_id": variation_id,
                                    "clinvar_allele_id": allele_id,
                                    "dbsnp_rs_id": rs_id,
                                    "gene_info": gene_info,
                                    "gene_symbols": gene_symbols,
                                    "gene_ids": gene_ids,
                                    "molecular_consequence": (
                                        molecular_consequence
                                    ),
                                    "clinical_hgvs": clinical_hgvs,
                                    "clinical_significance": (
                                        clinical_significance
                                    ),
                                    "label": label,
                                    "label_name": label_name,
                                    "review_status": review_status,
                                    "review_stars": review_stars,
                                    "variant_type": variant_type,
                                }
                            )

                            observed_counts["strict_records"] += 1
                            observed_counts[
                                f"{label_name}_records"
                            ] += 1
                            review_counts[
                                (review_stars, label_name)
                            ] += 1

        if len(genomic_keys) != observed_counts["strict_records"]:
            raise RuntimeError(
                "Unique genomic-key count does not equal exported "
                "record count."
            )

        verify_counts(observed_counts, review_counts)

        temporary_output.replace(OUTPUT_TSV)

    except Exception:
        if temporary_output.exists():
            temporary_output.unlink()

        raise

    sha256_checksum = calculate_sha256(OUTPUT_TSV)

    schema = build_schema(
        observed_counts,
        review_counts,
        annotation_missing_counts,
        sha256_checksum,
    )

    schema_text = (
        json.dumps(
            schema,
            indent=2,
            ensure_ascii=False,
        )
        + "\n"
    )

    checksum_text = (
        f"{sha256_checksum}  {OUTPUT_TSV.name}\n"
    )

    write_text_atomically(OUTPUT_SCHEMA, schema_text)
    write_text_atomically(OUTPUT_CHECKSUM, checksum_text)

    output_size = OUTPUT_TSV.stat().st_size

    print("CLINVAR BENCHMARK EXPORT")
    print("------------------------")
    print("Source VCF:", INPUT_VCF)
    print("Output TSV:", OUTPUT_TSV)
    print("Schema file:", OUTPUT_SCHEMA)
    print("Checksum file:", OUTPUT_CHECKSUM)
    print()

    print("EXPORTED RECORD COUNTS")
    print("----------------------")
    print("Total exported records:", observed_counts["strict_records"])
    print(
        "Pathogenic / likely pathogenic:",
        observed_counts["pathogenic_records"],
    )
    print(
        "Benign / likely benign:",
        observed_counts["benign_records"],
    )
    print("Unique genomic keys:", len(genomic_keys))
    print()

    print("REVIEW-STATUS COUNTS")
    print("--------------------")

    for stars in (2, 3, 4):
        pathogenic_count = review_counts[(stars, "pathogenic")]
        benign_count = review_counts[(stars, "benign")]

        print(
            f"{stars}-star pathogenic / likely pathogenic:",
            pathogenic_count,
        )
        print(
            f"{stars}-star benign / likely benign:",
            benign_count,
        )
        print(
            f"{stars}-star total:",
            pathogenic_count + benign_count,
        )

    print()
    print("OPTIONAL ANNOTATION MISSING COUNTS")
    print("----------------------------------")

    for annotation_name in (
        "clinvar_allele_id",
        "dbsnp_rs_id",
        "gene_info",
        "clinical_hgvs",
    ):
        print(
            f"{annotation_name}:",
            annotation_missing_counts[annotation_name],
        )

    print()
    print("FILE INTEGRITY")
    print("--------------")
    print("Compressed size in bytes:", output_size)
    print("SHA-256:", sha256_checksum)
    print()
    print(
        "PASS: Deterministic benchmark export completed for "
        "65,222 unique canonical strict records."
    )


if __name__ == "__main__":
    main()
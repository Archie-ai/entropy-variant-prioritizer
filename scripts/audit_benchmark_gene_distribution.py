import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


INPUT_FILE = Path(
    "data/processed/clinvar_benchmark_20260822.tsv.gz"
)

OUTPUT_DIRECTORY = Path("data/processed")

GENE_COUNTS_FILE = OUTPUT_DIRECTORY / (
    "clinvar_benchmark_20260822.gene_counts.tsv"
)

CHROMOSOME_COUNTS_FILE = OUTPUT_DIRECTORY / (
    "clinvar_benchmark_20260822.chromosome_counts.tsv"
)

MULTIGENE_RECORDS_FILE = OUTPUT_DIRECTORY / (
    "clinvar_benchmark_20260822.multigene_records.tsv"
)

SUMMARY_FILE = OUTPUT_DIRECTORY / (
    "clinvar_benchmark_20260822.gene_distribution_audit.json"
)

BENCHMARK_SHA256 = (
    "2a74db014bdf5830b1af7c4cd296c9ea1"
    "c8ea8d9e91859d1302987b265e7176f"
)

EXPECTED_RECORD_COUNT = 65222
EXPECTED_PATHOGENIC_COUNT = 24904
EXPECTED_BENIGN_COUNT = 40318

REQUIRED_COLUMNS = {
    "chrom",
    "pos",
    "ref",
    "alt",
    "genomic_key",
    "clinvar_variation_id",
    "gene_info",
    "gene_symbols",
    "gene_ids",
    "label",
    "label_name",
}


def chromosome_sort_key(chromosome):
    """Return a natural biological chromosome-sort key."""

    cleaned = chromosome.upper()

    if cleaned.startswith("CHR"):
        cleaned = cleaned[3:]

    if cleaned.isdigit():
        return 0, int(cleaned)

    special_chromosomes = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }

    if cleaned in special_chromosomes:
        return 0, special_chromosomes[cleaned]

    return 1, cleaned


def split_pipe_field(value):
    """Split a pipe-delimited field while preserving its order."""

    if not value:
        return []

    return [
        item.strip()
        for item in value.split("|")
    ]


def classify_gene_profile(pathogenic_count, benign_count):
    """Describe the clinical-label composition of one gene."""

    if pathogenic_count > 0 and benign_count > 0:
        return "mixed"

    if pathogenic_count > 0:
        return "pathogenic_only"

    if benign_count > 0:
        return "benign_only"

    return "no_valid_labels"


def gene_size_category(record_count):
    """Place a gene into a variant-count category."""

    if record_count == 1:
        return "1"

    if record_count <= 4:
        return "2-4"

    if record_count <= 9:
        return "5-9"

    if record_count <= 24:
        return "10-24"

    if record_count <= 49:
        return "25-49"

    if record_count <= 99:
        return "50-99"

    return "100+"


def write_tsv_atomically(
    output_path,
    fieldnames,
    rows,
):
    """Write a deterministic TSV through a temporary file."""

    temporary_path = output_path.with_name(
        output_path.name + ".tmp"
    )

    with temporary_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_MINIMAL,
        )

        writer.writeheader()
        writer.writerows(rows)

    temporary_path.replace(output_path)


def write_json_atomically(output_path, data):
    """Write deterministic JSON through a temporary file."""

    temporary_path = output_path.with_name(
        output_path.name + ".tmp"
    )

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

    temporary_path.replace(output_path)


def calculate_concentration(
    sorted_gene_counts,
    total_assignments,
    top_n,
):
    """Calculate the fraction of gene assignments in top genes."""

    assignments = sum(
        count
        for _, count in sorted_gene_counts[:top_n]
    )

    if total_assignments == 0:
        fraction = 0.0
    else:
        fraction = assignments / total_assignments

    return {
        "gene_count_requested": top_n,
        "gene_count_available": min(
            top_n,
            len(sorted_gene_counts),
        ),
        "assignments": assignments,
        "assignment_fraction": round(
            fraction,
            6,
        ),
        "assignment_percentage": round(
            fraction * 100,
            2,
        ),
    }


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Benchmark file was not found: {INPUT_FILE}"
        )

    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    total_records = 0
    class_counts = Counter()
    chromosome_counts = Counter()
    chromosome_class_counts = Counter()
    chromosome_multigene_counts = Counter()

    gene_total_counts = Counter()
    gene_class_counts = Counter()
    gene_chromosomes = defaultdict(set)

    symbol_to_gene_ids = defaultdict(set)
    gene_id_to_symbols = defaultdict(set)

    genomic_keys = set()
    variation_ids = set()

    single_gene_records = 0
    multigene_records = 0
    total_gene_assignments = 0

    records_with_missing_gene_id = 0
    records_with_repeated_gene_association = 0
    records_with_same_symbol_multiple_ids = 0
    records_with_same_id_multiple_symbols = 0

    multigene_output_rows = []

    with gzip.open(
        INPUT_FILE,
        "rt",
        encoding="utf-8",
        newline="",
    ) as file:
        reader = csv.DictReader(
            file,
            delimiter="\t",
        )

        observed_columns = set(
            reader.fieldnames or []
        )

        missing_columns = (
            REQUIRED_COLUMNS - observed_columns
        )

        if missing_columns:
            raise ValueError(
                "Benchmark is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for source_line_number, row in enumerate(
            reader,
            start=2,
        ):
            total_records += 1

            chrom = row["chrom"]
            genomic_key = row["genomic_key"]
            variation_id = row[
                "clinvar_variation_id"
            ]
            label_name = row["label_name"]

            if genomic_key in genomic_keys:
                raise ValueError(
                    "Duplicate genomic key at TSV line "
                    f"{source_line_number}: {genomic_key}"
                )

            genomic_keys.add(genomic_key)

            if variation_id in variation_ids:
                raise ValueError(
                    "Duplicate ClinVar Variation ID at TSV "
                    f"line {source_line_number}: "
                    f"{variation_id}"
                )

            variation_ids.add(variation_id)

            if label_name not in {
                "pathogenic",
                "benign",
            }:
                raise ValueError(
                    "Unexpected label name at TSV line "
                    f"{source_line_number}: "
                    f"{label_name!r}"
                )

            class_counts[label_name] += 1
            chromosome_counts[chrom] += 1
            chromosome_class_counts[
                (chrom, label_name)
            ] += 1

            gene_symbols = split_pipe_field(
                row["gene_symbols"]
            )
            gene_ids = split_pipe_field(
                row["gene_ids"]
            )

            if len(gene_symbols) != len(gene_ids):
                raise ValueError(
                    "Gene-symbol and Gene-ID counts differ "
                    f"at TSV line {source_line_number}."
                )

            associations = []

            for gene_symbol, gene_id in zip(
                gene_symbols,
                gene_ids,
            ):
                if not gene_symbol:
                    raise ValueError(
                        "Blank gene symbol at TSV line "
                        f"{source_line_number}."
                    )

                associations.append(
                    (gene_symbol, gene_id)
                )

            if not associations:
                raise ValueError(
                    "No usable gene association at TSV line "
                    f"{source_line_number}."
                )

            unique_associations = list(
                dict.fromkeys(associations)
            )

            if (
                len(unique_associations)
                != len(associations)
            ):
                records_with_repeated_gene_association += 1

            if any(
                not gene_id
                for _, gene_id in unique_associations
            ):
                records_with_missing_gene_id += 1

            row_symbol_to_ids = defaultdict(set)
            row_id_to_symbols = defaultdict(set)

            for gene_symbol, gene_id in (
                unique_associations
            ):
                row_symbol_to_ids[gene_symbol].add(
                    gene_id
                )

                if gene_id:
                    row_id_to_symbols[gene_id].add(
                        gene_symbol
                    )

            if any(
                len(gene_id_set) > 1
                for gene_id_set in (
                    row_symbol_to_ids.values()
                )
            ):
                records_with_same_symbol_multiple_ids += 1

            if any(
                len(symbol_set) > 1
                for symbol_set in (
                    row_id_to_symbols.values()
                )
            ):
                records_with_same_id_multiple_symbols += 1

            gene_count = len(unique_associations)

            if gene_count == 1:
                single_gene_records += 1
            else:
                multigene_records += 1
                chromosome_multigene_counts[
                    chrom
                ] += 1

                multigene_output_rows.append(
                    {
                        "genomic_key": genomic_key,
                        "clinvar_variation_id": (
                            variation_id
                        ),
                        "chrom": chrom,
                        "pos": row["pos"],
                        "ref": row["ref"],
                        "alt": row["alt"],
                        "label_name": label_name,
                        "gene_count": gene_count,
                        "gene_symbols": row[
                            "gene_symbols"
                        ],
                        "gene_ids": row["gene_ids"],
                    }
                )

            total_gene_assignments += gene_count

            for gene_symbol, gene_id in (
                unique_associations
            ):
                gene_key = (
                    gene_symbol,
                    gene_id,
                )

                gene_total_counts[gene_key] += 1
                gene_class_counts[
                    (
                        gene_symbol,
                        gene_id,
                        label_name,
                    )
                ] += 1

                gene_chromosomes[gene_key].add(
                    chrom
                )

                if gene_id:
                    symbol_to_gene_ids[
                        gene_symbol
                    ].add(gene_id)

                    gene_id_to_symbols[
                        gene_id
                    ].add(gene_symbol)

    if total_records != EXPECTED_RECORD_COUNT:
        raise RuntimeError(
            "Record count changed. Expected "
            f"{EXPECTED_RECORD_COUNT}; observed "
            f"{total_records}."
        )

    if (
        class_counts["pathogenic"]
        != EXPECTED_PATHOGENIC_COUNT
    ):
        raise RuntimeError(
            "Pathogenic count changed. Expected "
            f"{EXPECTED_PATHOGENIC_COUNT}; observed "
            f"{class_counts['pathogenic']}."
        )

    if (
        class_counts["benign"]
        != EXPECTED_BENIGN_COUNT
    ):
        raise RuntimeError(
            "Benign count changed. Expected "
            f"{EXPECTED_BENIGN_COUNT}; observed "
            f"{class_counts['benign']}."
        )

    if len(genomic_keys) != EXPECTED_RECORD_COUNT:
        raise RuntimeError(
            "Unique genomic-key count changed."
        )

    if len(variation_ids) != EXPECTED_RECORD_COUNT:
        raise RuntimeError(
            "Unique Variation-ID count changed."
        )

    if (
        single_gene_records + multigene_records
        != total_records
    ):
        raise RuntimeError(
            "Single-gene plus multigene records do "
            "not equal the benchmark total."
        )

    sorted_gene_counts = sorted(
        gene_total_counts.items(),
        key=lambda item: (
            -item[1],
            item[0][0],
            item[0][1],
        ),
    )

    gene_output_rows = []
    gene_profile_counts = Counter()
    gene_size_distribution = Counter()

    for gene_key, total_count in sorted_gene_counts:
        gene_symbol, gene_id = gene_key

        pathogenic_count = gene_class_counts[
            (
                gene_symbol,
                gene_id,
                "pathogenic",
            )
        ]

        benign_count = gene_class_counts[
            (
                gene_symbol,
                gene_id,
                "benign",
            )
        ]

        class_profile = classify_gene_profile(
            pathogenic_count,
            benign_count,
        )

        gene_profile_counts[class_profile] += 1

        gene_size_distribution[
            gene_size_category(total_count)
        ] += 1

        pathogenic_fraction = (
            pathogenic_count / total_count
        )

        benchmark_fraction = (
            total_count / total_records
        )

        chromosomes = sorted(
            gene_chromosomes[gene_key],
            key=chromosome_sort_key,
        )

        gene_output_rows.append(
            {
                "gene_symbol": gene_symbol,
                "gene_id": gene_id,
                "total_records": total_count,
                "pathogenic_records": (
                    pathogenic_count
                ),
                "benign_records": benign_count,
                "pathogenic_fraction": (
                    f"{pathogenic_fraction:.6f}"
                ),
                "benchmark_record_fraction": (
                    f"{benchmark_fraction:.6f}"
                ),
                "class_profile": class_profile,
                "chromosome_count": len(
                    chromosomes
                ),
                "chromosomes": "|".join(
                    chromosomes
                ),
            }
        )

    chromosome_output_rows = []

    sorted_chromosomes = sorted(
        chromosome_counts,
        key=chromosome_sort_key,
    )

    for chrom in sorted_chromosomes:
        total_count = chromosome_counts[chrom]
        pathogenic_count = (
            chromosome_class_counts[
                (chrom, "pathogenic")
            ]
        )
        benign_count = chromosome_class_counts[
            (chrom, "benign")
        ]

        chromosome_output_rows.append(
            {
                "chrom": chrom,
                "total_records": total_count,
                "pathogenic_records": (
                    pathogenic_count
                ),
                "benign_records": benign_count,
                "pathogenic_fraction": (
                    f"{pathogenic_count / total_count:.6f}"
                ),
                "benchmark_record_fraction": (
                    f"{total_count / total_records:.6f}"
                ),
                "multigene_records": (
                    chromosome_multigene_counts[
                        chrom
                    ]
                ),
            }
        )

    multigene_output_rows.sort(
        key=lambda row: (
            chromosome_sort_key(row["chrom"]),
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )
    )

    symbols_with_multiple_gene_ids = {
        gene_symbol: sorted(gene_ids)
        for gene_symbol, gene_ids in (
            symbol_to_gene_ids.items()
        )
        if len(gene_ids) > 1
    }

    gene_ids_with_multiple_symbols = {
        gene_id: sorted(gene_symbols)
        for gene_id, gene_symbols in (
            gene_id_to_symbols.items()
        )
        if len(gene_symbols) > 1
    }

    top_gene_rows = []

    for gene_key, total_count in (
        sorted_gene_counts[:25]
    ):
        gene_symbol, gene_id = gene_key

        top_gene_rows.append(
            {
                "gene_symbol": gene_symbol,
                "gene_id": gene_id,
                "total_records": total_count,
                "pathogenic_records": (
                    gene_class_counts[
                        (
                            gene_symbol,
                            gene_id,
                            "pathogenic",
                        )
                    ]
                ),
                "benign_records": (
                    gene_class_counts[
                        (
                            gene_symbol,
                            gene_id,
                            "benign",
                        )
                    ]
                ),
            }
        )

    summary = {
        "audit_name": (
            "ClinVar benchmark gene and chromosome "
            "distribution audit"
        ),
        "input_file": INPUT_FILE.as_posix(),
        "benchmark_sha256": BENCHMARK_SHA256,
        "record_count": total_records,
        "class_counts": {
            "pathogenic": (
                class_counts["pathogenic"]
            ),
            "benign": class_counts["benign"],
        },
        "gene_association_summary": {
            "single_gene_records": (
                single_gene_records
            ),
            "multigene_records": (
                multigene_records
            ),
            "total_gene_assignments": (
                total_gene_assignments
            ),
            "unique_gene_symbol_id_pairs": len(
                gene_total_counts
            ),
            "unique_gene_symbols": len(
                symbol_to_gene_ids
            ),
            "unique_nonmissing_gene_ids": len(
                gene_id_to_symbols
            ),
            "records_with_missing_gene_id": (
                records_with_missing_gene_id
            ),
            "records_with_repeated_gene_association": (
                records_with_repeated_gene_association
            ),
            "records_with_same_symbol_multiple_ids": (
                records_with_same_symbol_multiple_ids
            ),
            "records_with_same_id_multiple_symbols": (
                records_with_same_id_multiple_symbols
            ),
            "symbols_linked_to_multiple_gene_ids": (
                len(symbols_with_multiple_gene_ids)
            ),
            "gene_ids_linked_to_multiple_symbols": (
                len(gene_ids_with_multiple_symbols)
            ),
        },
        "gene_class_profiles": {
            "mixed": gene_profile_counts["mixed"],
            "pathogenic_only": gene_profile_counts[
                "pathogenic_only"
            ],
            "benign_only": gene_profile_counts[
                "benign_only"
            ],
        },
        "gene_record_count_distribution": {
            category: gene_size_distribution[
                category
            ]
            for category in (
                "1",
                "2-4",
                "5-9",
                "10-24",
                "25-49",
                "50-99",
                "100+",
            )
        },
        "gene_assignment_concentration": {
            "top_10": calculate_concentration(
                sorted_gene_counts,
                total_gene_assignments,
                10,
            ),
            "top_25": calculate_concentration(
                sorted_gene_counts,
                total_gene_assignments,
                25,
            ),
            "top_50": calculate_concentration(
                sorted_gene_counts,
                total_gene_assignments,
                50,
            ),
            "top_100": calculate_concentration(
                sorted_gene_counts,
                total_gene_assignments,
                100,
            ),
        },
        "chromosome_count": len(
            chromosome_counts
        ),
        "top_25_genes": top_gene_rows,
        "symbols_linked_to_multiple_gene_ids": (
            symbols_with_multiple_gene_ids
        ),
        "gene_ids_linked_to_multiple_symbols": (
            gene_ids_with_multiple_symbols
        ),
        "status": "PASS",
    }

    write_tsv_atomically(
        GENE_COUNTS_FILE,
        [
            "gene_symbol",
            "gene_id",
            "total_records",
            "pathogenic_records",
            "benign_records",
            "pathogenic_fraction",
            "benchmark_record_fraction",
            "class_profile",
            "chromosome_count",
            "chromosomes",
        ],
        gene_output_rows,
    )

    write_tsv_atomically(
        CHROMOSOME_COUNTS_FILE,
        [
            "chrom",
            "total_records",
            "pathogenic_records",
            "benign_records",
            "pathogenic_fraction",
            "benchmark_record_fraction",
            "multigene_records",
        ],
        chromosome_output_rows,
    )

    write_tsv_atomically(
        MULTIGENE_RECORDS_FILE,
        [
            "genomic_key",
            "clinvar_variation_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "label_name",
            "gene_count",
            "gene_symbols",
            "gene_ids",
        ],
        multigene_output_rows,
    )

    write_json_atomically(
        SUMMARY_FILE,
        summary,
    )

    print(
        "GENE AND CHROMOSOME DISTRIBUTION AUDIT"
    )
    print(
        "--------------------------------------"
    )
    print("Benchmark records:", total_records)
    print(
        "Pathogenic / likely pathogenic:",
        class_counts["pathogenic"],
    )
    print(
        "Benign / likely benign:",
        class_counts["benign"],
    )
    print()

    print("GENE-ASSOCIATION STRUCTURE")
    print("--------------------------")
    print(
        "Single-gene records:",
        single_gene_records,
    )
    print(
        "Multigene records:",
        multigene_records,
    )
    print(
        "Total gene assignments:",
        total_gene_assignments,
    )
    print(
        "Unique gene symbol/ID pairs:",
        len(gene_total_counts),
    )
    print(
        "Unique gene symbols:",
        len(symbol_to_gene_ids),
    )
    print(
        "Unique nonmissing Gene IDs:",
        len(gene_id_to_symbols),
    )
    print(
        "Records with missing Gene IDs:",
        records_with_missing_gene_id,
    )
    print()

    print("GENE CLASS PROFILES")
    print("-------------------")
    print(
        "Genes containing both classes:",
        gene_profile_counts["mixed"],
    )
    print(
        "Pathogenic-only genes:",
        gene_profile_counts[
            "pathogenic_only"
        ],
    )
    print(
        "Benign-only genes:",
        gene_profile_counts[
            "benign_only"
        ],
    )
    print()

    print("GENE-NAME MAPPING AUDIT")
    print("-----------------------")
    print(
        "Symbols linked to multiple Gene IDs:",
        len(symbols_with_multiple_gene_ids),
    )
    print(
        "Gene IDs linked to multiple symbols:",
        len(gene_ids_with_multiple_symbols),
    )
    print(
        "Records with repeated gene associations:",
        records_with_repeated_gene_association,
    )
    print()

    print("TOP 15 GENES BY RECORD COUNT")
    print("----------------------------")

    for rank, gene_row in enumerate(
        top_gene_rows[:15],
        start=1,
    ):
        print(
            f"{rank:>2}. "
            f"{gene_row['gene_symbol']} "
            f"(Gene ID {gene_row['gene_id']}): "
            f"{gene_row['total_records']} total; "
            f"{gene_row['pathogenic_records']} "
            f"pathogenic; "
            f"{gene_row['benign_records']} benign"
        )

    print()
    print("CHROMOSOME COUNTS")
    print("-----------------")

    for chromosome_row in chromosome_output_rows:
        print(
            f"{chromosome_row['chrom']}: "
            f"{chromosome_row['total_records']} total; "
            f"{chromosome_row['pathogenic_records']} "
            f"pathogenic; "
            f"{chromosome_row['benign_records']} benign; "
            f"{chromosome_row['multigene_records']} "
            f"multigene"
        )

    print()
    print("OUTPUT FILES")
    print("------------")
    print("Gene counts:", GENE_COUNTS_FILE)
    print(
        "Chromosome counts:",
        CHROMOSOME_COUNTS_FILE,
    )
    print(
        "Multigene records:",
        MULTIGENE_RECORDS_FILE,
    )
    print("JSON summary:", SUMMARY_FILE)
    print()

    print("FINAL AUDIT")
    print("-----------")
    print(
        "PASS: Gene and chromosome distributions "
        "were audited without changing the frozen "
        "65,222-record benchmark."
    )


if __name__ == "__main__":
    main()
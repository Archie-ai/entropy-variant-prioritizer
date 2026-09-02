import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


INPUT_FILE = Path(
    "data/processed/clinvar_benchmark_20260822.tsv.gz"
)

GENE_COMPONENTS_FILE = Path(
    "data/processed/"
    "clinvar_benchmark_20260822.gene_components.tsv"
)

GENE_ID_SYMBOL_FILE = Path(
    "data/processed/"
    "clinvar_benchmark_20260822.gene_id_symbol_audit.tsv"
)

JSON_SUMMARY_FILE = Path(
    "data/processed/"
    "clinvar_benchmark_20260822.gene_connectivity_audit.json"
)


EXPECTED_RECORDS = 65222
EXPECTED_PATHOGENIC = 24904
EXPECTED_BENIGN = 40318

EXPECTED_SINGLE_GENE_RECORDS = 58513
EXPECTED_MULTIGENE_RECORDS = 6709
EXPECTED_TOTAL_GENE_ASSIGNMENTS = 72408

EXPECTED_GENE_SYMBOL_ID_PAIRS = 8996
EXPECTED_UNIQUE_GENE_SYMBOLS = 8996
EXPECTED_UNIQUE_GENE_IDS = 8989

EXPECTED_GENE_IDS_WITH_MULTIPLE_SYMBOLS = 7
EXPECTED_SYMBOLS_WITH_MULTIPLE_GENE_IDS = 0
EXPECTED_REPEATED_GENE_ID_RECORDS = 0


class UnionFind:
    """Maintain connected groups of Gene IDs."""

    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, item):
        if item not in self.parent:
            self.parent[item] = item
            self.rank[item] = 0

    def find(self, item):
        parent = self.parent[item]

        if parent != item:
            self.parent[item] = self.find(parent)

        return self.parent[item]

    def union(self, item_a, item_b):
        root_a = self.find(item_a)
        root_b = self.find(item_b)

        if root_a == root_b:
            return

        rank_a = self.rank[root_a]
        rank_b = self.rank[root_b]

        if rank_a < rank_b:
            self.parent[root_a] = root_b

        elif rank_a > rank_b:
            self.parent[root_b] = root_a

        else:
            self.parent[root_b] = root_a
            self.rank[root_a] += 1


def gene_id_sort_key(gene_id):
    """Sort numeric Gene IDs numerically and other IDs alphabetically."""

    if gene_id.isdigit():
        return 0, int(gene_id)

    return 1, gene_id


def parse_label(raw_label):
    """Convert the benchmark label into a readable class name."""

    value = raw_label.strip().lower()

    pathogenic_values = {
        "1",
        "pathogenic",
        "pathogenic_likely_pathogenic",
        "pathogenic / likely pathogenic",
    }

    benign_values = {
        "0",
        "benign",
        "benign_likely_benign",
        "benign / likely benign",
    }

    if value in pathogenic_values:
        return "pathogenic"

    if value in benign_values:
        return "benign"

    raise ValueError(
        f"Unrecognized benchmark label: {raw_label!r}"
    )


def parse_gene_info(raw_gene_info):
    """
    Parse ClinVar GENEINFO values.

    Expected format:
        SYMBOL:GENE_ID

    Multigene example:
        SYMBOL1:GENE_ID1|SYMBOL2:GENE_ID2
    """

    value = raw_gene_info.strip()

    if not value or value in {".", "-"}:
        raise ValueError(
            "A benchmark record has missing gene_info."
        )

    associations = []

    for entry in value.split("|"):
        entry = entry.strip()

        if not entry:
            continue

        if ":" not in entry:
            raise ValueError(
                f"Malformed gene_info entry: {entry!r}"
            )

        symbol, gene_id = entry.rsplit(":", 1)

        symbol = symbol.strip()
        gene_id = gene_id.strip()

        if not symbol:
            raise ValueError(
                f"Missing gene symbol in entry: {entry!r}"
            )

        if not gene_id:
            raise ValueError(
                f"Missing Gene ID in entry: {entry!r}"
            )

        associations.append((symbol, gene_id))

    if not associations:
        raise ValueError(
            f"No valid associations in gene_info: "
            f"{raw_gene_info!r}"
        )

    return associations


def component_profile(pathogenic_count, benign_count):
    """Describe which benchmark classes occur in a component."""

    if pathogenic_count > 0 and benign_count > 0:
        return "both_classes"

    if pathogenic_count > 0:
        return "pathogenic_only"

    if benign_count > 0:
        return "benign_only"

    return "no_records"


def component_size_bin(gene_count):
    """Assign a connected component to a size category."""

    if gene_count == 1:
        return "1"

    if gene_count == 2:
        return "2"

    if 3 <= gene_count <= 5:
        return "3-5"

    if 6 <= gene_count <= 10:
        return "6-10"

    if 11 <= gene_count <= 50:
        return "11-50"

    return "51+"


def shorten_list(values, limit=12):
    """Create a compact representation for terminal output."""

    values = list(values)

    if len(values) <= limit:
        return "|".join(values)

    visible = "|".join(values[:limit])
    hidden_count = len(values) - limit

    return f"{visible}|...(+{hidden_count})"


if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Input benchmark was not found: {INPUT_FILE}"
    )


records = []

record_count = 0
pathogenic_count = 0
benign_count = 0

single_gene_records = 0
multigene_records = 0
total_gene_assignments = 0
repeated_gene_id_records = 0

gene_symbol_id_pairs = set()
gene_symbols = set()
gene_ids = set()

gene_id_to_symbols = defaultdict(set)
symbol_to_gene_ids = defaultdict(set)

gene_record_counts = Counter()
gene_pathogenic_counts = Counter()
gene_benign_counts = Counter()

union_find = UnionFind()


with gzip.open(
    INPUT_FILE,
    "rt",
    encoding="utf-8",
    newline="",
) as input_handle:

    reader = csv.DictReader(
        input_handle,
        delimiter="\t",
    )

    if reader.fieldnames is None:
        raise ValueError(
            "The benchmark file does not contain a header."
        )

    required_columns = {
        "label",
        "gene_info",
    }

    missing_columns = (
        required_columns - set(reader.fieldnames)
    )

    if missing_columns:
        raise ValueError(
            "Missing required benchmark columns: "
            + ", ".join(sorted(missing_columns))
        )

    for row_number, row in enumerate(reader, start=2):
        label = parse_label(row["label"])

        try:
            associations = parse_gene_info(
                row["gene_info"]
            )

        except ValueError as error:
            raise ValueError(
                f"Row {row_number}: {error}"
            ) from error

        record_count += 1
        total_gene_assignments += len(associations)

        if label == "pathogenic":
            pathogenic_count += 1

        else:
            benign_count += 1

        record_gene_ids = []

        for symbol, gene_id in associations:
            gene_symbol_id_pairs.add(
                (symbol, gene_id)
            )

            gene_symbols.add(symbol)
            gene_ids.add(gene_id)

            gene_id_to_symbols[gene_id].add(symbol)
            symbol_to_gene_ids[symbol].add(gene_id)

            union_find.add(gene_id)
            record_gene_ids.append(gene_id)

        unique_record_gene_ids = sorted(
            set(record_gene_ids),
            key=gene_id_sort_key,
        )

        if len(unique_record_gene_ids) != len(
            record_gene_ids
        ):
            repeated_gene_id_records += 1

        if len(unique_record_gene_ids) == 1:
            single_gene_records += 1

        else:
            multigene_records += 1

            first_gene_id = unique_record_gene_ids[0]

            for other_gene_id in unique_record_gene_ids[1:]:
                union_find.union(
                    first_gene_id,
                    other_gene_id,
                )

        for gene_id in unique_record_gene_ids:
            gene_record_counts[gene_id] += 1

            if label == "pathogenic":
                gene_pathogenic_counts[gene_id] += 1

            else:
                gene_benign_counts[gene_id] += 1

        records.append(
            {
                "label": label,
                "gene_ids": unique_record_gene_ids,
            }
        )


gene_ids_with_multiple_symbols = {
    gene_id: sorted(symbols)
    for gene_id, symbols in gene_id_to_symbols.items()
    if len(symbols) > 1
}

symbols_with_multiple_gene_ids = {
    symbol: sorted(
        linked_gene_ids,
        key=gene_id_sort_key,
    )
    for symbol, linked_gene_ids
    in symbol_to_gene_ids.items()
    if len(linked_gene_ids) > 1
}


root_to_gene_ids = defaultdict(list)

for gene_id in sorted(
    gene_ids,
    key=gene_id_sort_key,
):
    root = union_find.find(gene_id)
    root_to_gene_ids[root].append(gene_id)


component_members = {}
gene_to_component = {}

for member_gene_ids in root_to_gene_ids.values():
    sorted_gene_ids = sorted(
        member_gene_ids,
        key=gene_id_sort_key,
    )

    minimum_gene_id = sorted_gene_ids[0]
    component_id = f"GC_{minimum_gene_id}"

    component_members[component_id] = sorted_gene_ids

    for gene_id in sorted_gene_ids:
        gene_to_component[gene_id] = component_id


component_record_counts = Counter()
component_pathogenic_counts = Counter()
component_benign_counts = Counter()

records_with_invalid_component_mapping = 0

for record in records:
    record_component_ids = {
        gene_to_component[gene_id]
        for gene_id in record["gene_ids"]
    }

    if len(record_component_ids) != 1:
        records_with_invalid_component_mapping += 1
        continue

    component_id = next(iter(record_component_ids))

    component_record_counts[component_id] += 1

    if record["label"] == "pathogenic":
        component_pathogenic_counts[component_id] += 1

    else:
        component_benign_counts[component_id] += 1


component_rows = []

for component_id, member_gene_ids in (
    component_members.items()
):
    component_symbols = sorted(
        {
            symbol
            for gene_id in member_gene_ids
            for symbol in gene_id_to_symbols[gene_id]
        }
    )

    component_pathogenic = (
        component_pathogenic_counts[component_id]
    )

    component_benign = (
        component_benign_counts[component_id]
    )

    component_rows.append(
        {
            "gene_component_id": component_id,
            "gene_count": len(member_gene_ids),
            "record_count":
                component_record_counts[component_id],
            "pathogenic_count":
                component_pathogenic,
            "benign_count":
                component_benign,
            "class_profile": component_profile(
                component_pathogenic,
                component_benign,
            ),
            "gene_ids": "|".join(member_gene_ids),
            "gene_symbols": "|".join(
                component_symbols
            ),
        }
    )


component_rows.sort(
    key=lambda row: (
        -row["record_count"],
        -row["gene_count"],
        row["gene_component_id"],
    )
)


single_gene_components = sum(
    row["gene_count"] == 1
    for row in component_rows
)

multigene_components = sum(
    row["gene_count"] > 1
    for row in component_rows
)

largest_component_gene_count = max(
    row["gene_count"]
    for row in component_rows
)

largest_component_record_count = max(
    row["record_count"]
    for row in component_rows
)

records_in_single_gene_components = sum(
    row["record_count"]
    for row in component_rows
    if row["gene_count"] == 1
)

records_in_multigene_components = sum(
    row["record_count"]
    for row in component_rows
    if row["gene_count"] > 1
)

components_with_both_classes = sum(
    row["class_profile"] == "both_classes"
    for row in component_rows
)

pathogenic_only_components = sum(
    row["class_profile"] == "pathogenic_only"
    for row in component_rows
)

benign_only_components = sum(
    row["class_profile"] == "benign_only"
    for row in component_rows
)


size_bin_order = [
    "1",
    "2",
    "3-5",
    "6-10",
    "11-50",
    "51+",
]

component_size_distribution = {
    size_bin: {
        "component_count": 0,
        "gene_count": 0,
        "record_count": 0,
        "pathogenic_count": 0,
        "benign_count": 0,
    }
    for size_bin in size_bin_order
}

for row in component_rows:
    size_bin = component_size_bin(
        row["gene_count"]
    )

    values = component_size_distribution[size_bin]

    values["component_count"] += 1
    values["gene_count"] += row["gene_count"]
    values["record_count"] += row["record_count"]
    values["pathogenic_count"] += (
        row["pathogenic_count"]
    )
    values["benign_count"] += row["benign_count"]


assert record_count == EXPECTED_RECORDS
assert pathogenic_count == EXPECTED_PATHOGENIC
assert benign_count == EXPECTED_BENIGN

assert (
    single_gene_records
    == EXPECTED_SINGLE_GENE_RECORDS
)

assert (
    multigene_records
    == EXPECTED_MULTIGENE_RECORDS
)

assert (
    total_gene_assignments
    == EXPECTED_TOTAL_GENE_ASSIGNMENTS
)

assert (
    len(gene_symbol_id_pairs)
    == EXPECTED_GENE_SYMBOL_ID_PAIRS
)

assert (
    len(gene_symbols)
    == EXPECTED_UNIQUE_GENE_SYMBOLS
)

assert (
    len(gene_ids)
    == EXPECTED_UNIQUE_GENE_IDS
)

assert (
    len(gene_ids_with_multiple_symbols)
    == EXPECTED_GENE_IDS_WITH_MULTIPLE_SYMBOLS
)

assert (
    len(symbols_with_multiple_gene_ids)
    == EXPECTED_SYMBOLS_WITH_MULTIPLE_GENE_IDS
)

assert (
    repeated_gene_id_records
    == EXPECTED_REPEATED_GENE_ID_RECORDS
)

assert records_with_invalid_component_mapping == 0

assert sum(
    row["record_count"]
    for row in component_rows
) == EXPECTED_RECORDS

assert sum(
    row["pathogenic_count"]
    for row in component_rows
) == EXPECTED_PATHOGENIC

assert sum(
    row["benign_count"]
    for row in component_rows
) == EXPECTED_BENIGN


GENE_COMPONENTS_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)


with GENE_COMPONENTS_FILE.open(
    "w",
    encoding="utf-8",
    newline="",
) as output_handle:

    fieldnames = [
        "gene_component_id",
        "gene_count",
        "record_count",
        "pathogenic_count",
        "benign_count",
        "class_profile",
        "gene_ids",
        "gene_symbols",
    ]

    writer = csv.DictWriter(
        output_handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(component_rows)


with GENE_ID_SYMBOL_FILE.open(
    "w",
    encoding="utf-8",
    newline="",
) as output_handle:

    fieldnames = [
        "gene_id",
        "symbol_count",
        "symbols",
        "record_count",
        "pathogenic_count",
        "benign_count",
        "gene_component_id",
        "component_gene_count",
    ]

    writer = csv.DictWriter(
        output_handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for gene_id in sorted(
        gene_ids,
        key=gene_id_sort_key,
    ):
        component_id = gene_to_component[gene_id]

        writer.writerow(
            {
                "gene_id": gene_id,
                "symbol_count": len(
                    gene_id_to_symbols[gene_id]
                ),
                "symbols": "|".join(
                    sorted(
                        gene_id_to_symbols[gene_id]
                    )
                ),
                "record_count":
                    gene_record_counts[gene_id],
                "pathogenic_count":
                    gene_pathogenic_counts[gene_id],
                "benign_count":
                    gene_benign_counts[gene_id],
                "gene_component_id":
                    component_id,
                "component_gene_count": len(
                    component_members[component_id]
                ),
            }
        )


top_components_for_json = []

for row in component_rows[:25]:
    top_components_for_json.append(
        {
            "gene_component_id":
                row["gene_component_id"],
            "gene_count": row["gene_count"],
            "record_count": row["record_count"],
            "pathogenic_count":
                row["pathogenic_count"],
            "benign_count":
                row["benign_count"],
            "class_profile":
                row["class_profile"],
            "gene_ids":
                row["gene_ids"].split("|"),
            "gene_symbols":
                row["gene_symbols"].split("|"),
        }
    )


json_summary = {
    "input_file": str(INPUT_FILE),
    "benchmark": {
        "records": record_count,
        "pathogenic": pathogenic_count,
        "benign": benign_count,
    },
    "gene_associations": {
        "single_gene_records":
            single_gene_records,
        "multigene_records":
            multigene_records,
        "total_gene_assignments":
            total_gene_assignments,
        "unique_gene_symbol_id_pairs":
            len(gene_symbol_id_pairs),
        "unique_gene_symbols":
            len(gene_symbols),
        "unique_gene_ids":
            len(gene_ids),
    },
    "gene_name_mapping": {
        "gene_ids_with_multiple_symbols":
            gene_ids_with_multiple_symbols,
        "symbols_with_multiple_gene_ids":
            symbols_with_multiple_gene_ids,
        "repeated_gene_id_records":
            repeated_gene_id_records,
    },
    "connectivity": {
        "total_components":
            len(component_rows),
        "single_gene_components":
            single_gene_components,
        "multigene_components":
            multigene_components,
        "largest_component_gene_count":
            largest_component_gene_count,
        "largest_component_record_count":
            largest_component_record_count,
        "records_in_single_gene_components":
            records_in_single_gene_components,
        "records_in_multigene_components":
            records_in_multigene_components,
        "components_with_both_classes":
            components_with_both_classes,
        "pathogenic_only_components":
            pathogenic_only_components,
        "benign_only_components":
            benign_only_components,
        "invalid_component_mappings":
            records_with_invalid_component_mapping,
    },
    "component_size_distribution":
        component_size_distribution,
    "top_components_by_record_count":
        top_components_for_json,
    "output_files": {
        "gene_components":
            str(GENE_COMPONENTS_FILE),
        "gene_id_symbol_audit":
            str(GENE_ID_SYMBOL_FILE),
        "json_summary":
            str(JSON_SUMMARY_FILE),
    },
}


with JSON_SUMMARY_FILE.open(
    "w",
    encoding="utf-8",
) as output_handle:

    json.dump(
        json_summary,
        output_handle,
        indent=2,
        sort_keys=True,
    )

    output_handle.write("\n")


multigene_percentage = (
    multigene_records / record_count * 100
)

average_genes_per_multigene_record = (
    (
        total_gene_assignments
        - single_gene_records
    )
    / multigene_records
)


print("GENE CONNECTIVITY AUDIT")
print("-----------------------")
print("Benchmark records:", record_count)
print(
    "Pathogenic / likely pathogenic:",
    pathogenic_count,
)
print(
    "Benign / likely benign:",
    benign_count,
)

print()
print("MULTIGENE RECORD SUMMARY")
print("------------------------")
print("Single-gene records:", single_gene_records)
print("Multigene records:", multigene_records)
print(
    "Multigene percentage:",
    f"{multigene_percentage:.2f}%",
)
print(
    "Average genes per multigene record:",
    f"{average_genes_per_multigene_record:.2f}",
)
print(
    "Total gene assignments:",
    total_gene_assignments,
)

print()
print("GENE ID TO SYMBOL AUDIT")
print("-----------------------")
print(
    "Unique Gene IDs:",
    len(gene_ids),
)
print(
    "Unique gene symbols:",
    len(gene_symbols),
)
print(
    "Gene IDs linked to multiple symbols:",
    len(gene_ids_with_multiple_symbols),
)
print(
    "Symbols linked to multiple Gene IDs:",
    len(symbols_with_multiple_gene_ids),
)

print()
print("GENE IDs LINKED TO MULTIPLE SYMBOLS")
print("-----------------------------------")

if gene_ids_with_multiple_symbols:
    for gene_id in sorted(
        gene_ids_with_multiple_symbols,
        key=gene_id_sort_key,
    ):
        symbols = gene_ids_with_multiple_symbols[
            gene_id
        ]

        print(
            f"Gene ID {gene_id}: "
            + " | ".join(symbols)
        )

else:
    print("None")

print()
print("CONNECTED-COMPONENT STRUCTURE")
print("-----------------------------")
print(
    "Total gene components:",
    len(component_rows),
)
print(
    "Single-gene components:",
    single_gene_components,
)
print(
    "Multigene connected components:",
    multigene_components,
)
print(
    "Largest component, genes:",
    largest_component_gene_count,
)
print(
    "Largest component, records:",
    largest_component_record_count,
)
print(
    "Records in single-gene components:",
    records_in_single_gene_components,
)
print(
    "Records in multigene components:",
    records_in_multigene_components,
)

print()
print("COMPONENT CLASS PROFILES")
print("------------------------")
print(
    "Components containing both classes:",
    components_with_both_classes,
)
print(
    "Pathogenic-only components:",
    pathogenic_only_components,
)
print(
    "Benign-only components:",
    benign_only_components,
)

print()
print("COMPONENT SIZE DISTRIBUTION")
print("---------------------------")

for size_bin in size_bin_order:
    values = component_size_distribution[
        size_bin
    ]

    print(
        f"{size_bin} genes: "
        f"{values['component_count']} components; "
        f"{values['gene_count']} genes; "
        f"{values['record_count']} records; "
        f"{values['pathogenic_count']} pathogenic; "
        f"{values['benign_count']} benign"
    )

print()
print("TOP 15 COMPONENTS BY RECORD COUNT")
print("---------------------------------")

for rank, row in enumerate(
    component_rows[:15],
    start=1,
):
    symbols = row["gene_symbols"].split("|")

    print(
        f"{rank:2d}. "
        f"{row['gene_component_id']}: "
        f"{row['gene_count']} genes; "
        f"{row['record_count']} records; "
        f"{row['pathogenic_count']} pathogenic; "
        f"{row['benign_count']} benign; "
        f"{shorten_list(symbols)}"
    )

print()
print("OUTPUT FILES")
print("------------")
print(
    "Gene components:",
    GENE_COMPONENTS_FILE,
)
print(
    "Gene ID-symbol audit:",
    GENE_ID_SYMBOL_FILE,
)
print(
    "JSON summary:",
    JSON_SUMMARY_FILE,
)

print()
print("FINAL AUDIT")
print("-----------")
print(
    "PASS: Every benchmark record maps to exactly one "
    "gene-connected component, without changing the "
    "frozen 65,222-record benchmark."
)
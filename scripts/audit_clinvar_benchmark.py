import gzip
from collections import Counter, defaultdict


VCF_FILE = "data/raw/clinvar_20260822.vcf.gz"
EXPECTED_STRICT_COUNT = 65222

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

ACCEPTED_REVIEW_STATUSES = {
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
}

REVIEW_STATUS_LABELS = {
    "criteria_provided,_multiple_submitters,_no_conflicts": "2-star",
    "reviewed_by_expert_panel": "3-star",
    "practice_guideline": "4-star",
}


# Benchmark reproduction counters
total_variants = 0
missense_variants = 0
snv_variants = 0
one_character_variants = 0
noncanonical_variants = 0
canonical_variants = 0
clinical_significance_variants = 0
strict_variants = 0

pathogenic_variants = 0
benign_variants = 0

strict_pathogenic_variants = 0
strict_benign_variants = 0

# VCF and allele-integrity counters
malformed_records = 0
multiallelic_records = 0
identical_ref_alt_records = 0
strict_records_missing_clinvar_id = 0

# Duplicate and conflict tracking
genomic_key_counts = Counter()
genomic_key_clinvar_ids = defaultdict(set)
genomic_key_labels = defaultdict(set)
clinvar_id_counts = Counter()

# Review-status tracking
review_status_class_counts = defaultdict(Counter)


with gzip.open(VCF_FILE, "rt") as file:
    for line in file:

        # Skip metadata and the VCF column-header line
        if line.startswith("#"):
            continue

        total_variants += 1

        fields = line.rstrip("\n").split("\t")

        # A valid VCF row must have at least eight fixed columns
        if len(fields) < 8:
            malformed_records += 1
            continue

        chrom = fields[0]
        position = fields[1]
        clinvar_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        info = fields[7]

        # Parse the INFO field into a dictionary
        info_dict = {}

        for item in info.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info_dict[key] = value
            else:
                info_dict[item] = True

        molecular_consequence = info_dict.get("MC", "")
        variant_type = info_dict.get("CLNVC", "")
        clinical_significance = info_dict.get("CLNSIG", "")
        review_status = info_dict.get("CLNREVSTAT", "")

        # Filter 1: missense consequence
        if "missense_variant" not in molecular_consequence:
            continue

        missense_variants += 1

        # Filter 2: ClinVar SNV annotation
        if variant_type != "single_nucleotide_variant":
            continue

        snv_variants += 1

        # Inspect multiallelic ALT values among annotated missense SNVs
        if "," in alt:
            multiallelic_records += 1

        # Filter 3: REF and ALT must each contain exactly one character
        if len(ref) != 1 or len(alt) != 1:
            continue

        one_character_variants += 1

        # Inspect records in which REF and ALT are identical
        if ref == alt:
            identical_ref_alt_records += 1

        # Filter 4: REF and ALT must be canonical DNA nucleotides
        if ref not in CANONICAL_BASES or alt not in CANONICAL_BASES:
            noncanonical_variants += 1
            continue

        canonical_variants += 1

        # Filter 5: accepted binary clinical significance
        if clinical_significance in PATHOGENIC_CLASSES:
            label = "pathogenic"
            pathogenic_variants += 1

        elif clinical_significance in BENIGN_CLASSES:
            label = "benign"
            benign_variants += 1

        else:
            continue

        clinical_significance_variants += 1

        # Filter 6: strict 2-4-star review status
        if review_status not in ACCEPTED_REVIEW_STATUSES:
            continue

        strict_variants += 1

        if label == "pathogenic":
            strict_pathogenic_variants += 1
        else:
            strict_benign_variants += 1

        # Record exact review-status and class counts
        review_status_class_counts[review_status][label] += 1

        # Construct the genomic variant key
        genomic_key = (chrom, position, ref, alt)

        genomic_key_counts[genomic_key] += 1
        genomic_key_labels[genomic_key].add(label)

        # Audit ClinVar Variation IDs
        if clinvar_id in {"", "."}:
            strict_records_missing_clinvar_id += 1
        else:
            genomic_key_clinvar_ids[genomic_key].add(clinvar_id)
            clinvar_id_counts[clinvar_id] += 1


# Calculate duplicate statistics
unique_genomic_keys = len(genomic_key_counts)

duplicate_genomic_keys = sum(
    1
    for count in genomic_key_counts.values()
    if count > 1
)

duplicate_record_excess = sum(
    count - 1
    for count in genomic_key_counts.values()
    if count > 1
)

maximum_records_per_key = max(
    genomic_key_counts.values(),
    default=0,
)

keys_with_multiple_clinvar_ids = sum(
    1
    for clinvar_ids in genomic_key_clinvar_ids.values()
    if len(clinvar_ids) > 1
)

clinvar_ids_in_multiple_records = sum(
    1
    for count in clinvar_id_counts.values()
    if count > 1
)

# Identify cross-label conflicts
cross_label_conflict_keys = [
    genomic_key
    for genomic_key, labels in genomic_key_labels.items()
    if len(labels) > 1
]

# Identify duplicate examples
duplicate_examples = sorted(
    (
        (genomic_key, count)
        for genomic_key, count in genomic_key_counts.items()
        if count > 1
    ),
    key=lambda item: (-item[1], item[0]),
)[:10]


print("BENCHMARK COUNT REPRODUCTION")
print("----------------------------")
print("Total ClinVar variants:", total_variants)
print("Missense variants:", missense_variants)
print("Missense SNVs:", snv_variants)
print("One-character REF/ALT SNVs:", one_character_variants)
print("Excluded noncanonical one-character alleles:", noncanonical_variants)
print("Canonical A/C/G/T missense SNVs:", canonical_variants)
print("Accepted clinical significance:", clinical_significance_variants)
print("Strict 2-4-star records:", strict_variants)
print(
    "Strict pathogenic / likely pathogenic:",
    strict_pathogenic_variants,
)
print(
    "Strict benign / likely benign:",
    strict_benign_variants,
)
print(
    "Class-count verification:",
    strict_pathogenic_variants + strict_benign_variants,
)


print("\nVCF AND ALLELE INTEGRITY")
print("------------------------")
print("Malformed VCF records:", malformed_records)
print(
    "Multiallelic ALT records among ClinVar-annotated missense SNVs:",
    multiallelic_records,
)
print(
    "Noncanonical one-character REF/ALT alleles excluded:",
    noncanonical_variants,
)
print(
    "Records with identical REF and ALT:",
    identical_ref_alt_records,
)
print(
    "Strict records missing a ClinVar Variation ID:",
    strict_records_missing_clinvar_id,
)


print("\nUNIQUE VARIANT AND DUPLICATE AUDIT")
print("----------------------------------")
print("Strict candidate records:", strict_variants)
print("Unique CHROM/POS/REF/ALT keys:", unique_genomic_keys)
print(
    "Genomic keys occurring more than once:",
    duplicate_genomic_keys,
)
print("Duplicate-record excess:", duplicate_record_excess)
print(
    "Unique-key plus duplicate-excess verification:",
    unique_genomic_keys + duplicate_record_excess,
)
print(
    "Maximum records for one genomic key:",
    maximum_records_per_key,
)
print(
    "Keys linked to multiple ClinVar IDs:",
    keys_with_multiple_clinvar_ids,
)
print(
    "ClinVar IDs appearing in multiple strict records:",
    clinvar_ids_in_multiple_records,
)


print("\nCROSS-LABEL CONFLICT AUDIT")
print("--------------------------")
print(
    "Genomic keys with both pathogenic and benign labels:",
    len(cross_label_conflict_keys),
)


print("\nEXACT REVIEW-STATUS BREAKDOWN")
print("-----------------------------")

review_status_order = [
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
]

review_status_total = 0

for review_status in review_status_order:
    status_label = REVIEW_STATUS_LABELS[review_status]

    pathogenic_count = review_status_class_counts[
        review_status
    ]["pathogenic"]

    benign_count = review_status_class_counts[
        review_status
    ]["benign"]

    status_total = pathogenic_count + benign_count
    review_status_total += status_total

    if strict_variants > 0:
        percentage = (status_total / strict_variants) * 100
    else:
        percentage = 0.0

    print(f"\n{status_label}")
    print("Review status:", review_status)
    print(
        "Pathogenic / likely pathogenic:",
        pathogenic_count,
    )
    print(
        "Benign / likely benign:",
        benign_count,
    )
    print("Total:", status_total)
    print(
        "Percentage of strict benchmark:",
        f"{percentage:.2f}%",
    )


print("\nREVIEW-STATUS VERIFICATION")
print("--------------------------")
print(
    "Sum of 2-star, 3-star, and 4-star records:",
    review_status_total,
)


print("\nTOP DUPLICATE GENOMIC-KEY EXAMPLES")
print("----------------------------------")

if duplicate_examples:
    for genomic_key, count in duplicate_examples:
        chrom, position, ref, alt = genomic_key

        print(
            chrom,
            position,
            ref,
            alt,
            count,
            sep="\t",
        )
else:
    print("No duplicate genomic keys were found.")


print("\nCROSS-LABEL CONFLICT EXAMPLES")
print("-----------------------------")

if cross_label_conflict_keys:
    for genomic_key in cross_label_conflict_keys[:10]:
        chrom, position, ref, alt = genomic_key
        labels = sorted(genomic_key_labels[genomic_key])

        print(
            chrom,
            position,
            ref,
            alt,
            ",".join(labels),
            sep="\t",
        )
else:
    print("No cross-label conflicts were found.")


# Final accounting and integrity checks
accounting_checks = [
    strict_variants == EXPECTED_STRICT_COUNT,
    strict_variants
    == strict_pathogenic_variants + strict_benign_variants,
    clinical_significance_variants
    == pathogenic_variants + benign_variants,
    one_character_variants
    == canonical_variants + noncanonical_variants,
    strict_variants
    == unique_genomic_keys + duplicate_record_excess,
    strict_variants == review_status_total,
]

integrity_checks = [
    malformed_records == 0,
    multiallelic_records == 0,
    identical_ref_alt_records == 0,
    strict_records_missing_clinvar_id == 0,
    duplicate_genomic_keys == 0,
    keys_with_multiple_clinvar_ids == 0,
    clinvar_ids_in_multiple_records == 0,
    len(cross_label_conflict_keys) == 0,
]


print("\nFINAL INTEGRITY CHECK")
print("---------------------")

if all(accounting_checks) and all(integrity_checks):
    print(
        "PASS: All accounting and integrity checks passed "
        f"for {strict_variants:,} canonical strict records."
    )
else:
    print("FAIL: One or more integrity checks did not pass.")

    if strict_variants != EXPECTED_STRICT_COUNT:
        print(
            "Expected strict count:",
            EXPECTED_STRICT_COUNT,
        )
        print(
            "Observed strict count:",
            strict_variants,
        )
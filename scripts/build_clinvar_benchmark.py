import gzip


vcf_file = "data/raw/clinvar_20260822.vcf.gz"

# Filtering-stage counters
total_variants = 0
missense_variants = 0
snv_variants = 0
allele_valid_variants = 0
clinical_significance_variants = 0
high_review_variants = 0

# Class counters before the review-status filter
pathogenic_variants = 0
benign_variants = 0

# Class counters after the 2-4-star review-status filter
high_review_pathogenic_variants = 0
high_review_benign_variants = 0


pathogenic_classes = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic"
}

benign_classes = {
    "Benign",
    "Likely_benign",
    "Benign/Likely_benign"
}

# ClinVar review statuses corresponding to 2, 3, or 4 stars
accepted_review_statuses = {
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline"
}


with gzip.open(vcf_file, "rt") as file:

    for line in file:

        # Skip VCF metadata and column-header lines
        if line.startswith("#"):
            continue

        total_variants += 1

        # Split the VCF row into columns
        fields = line.strip().split("\t")

        ref = fields[3]
        alt = fields[4]
        info = fields[7]

        # Parse the INFO column into a dictionary
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

        # Filter 1: retain missense variants
        if "missense_variant" not in molecular_consequence:
            continue

        missense_variants += 1

        # Filter 2: retain variants annotated by ClinVar as SNVs
        if variant_type != "single_nucleotide_variant":
            continue

        snv_variants += 1

        # Filter 3: REF and ALT must each contain one nucleotide
        if len(ref) != 1 or len(alt) != 1:
            continue

        allele_valid_variants += 1

        # Filter 4: assign accepted binary clinical labels
        if clinical_significance in pathogenic_classes:
            label = 1
            pathogenic_variants += 1

        elif clinical_significance in benign_classes:
            label = 0
            benign_variants += 1

        else:
            continue

        clinical_significance_variants += 1

        # Filter 5: retain ClinVar records with 2-4-star review status
        if review_status not in accepted_review_statuses:
            continue

        high_review_variants += 1

        # Count pathogenic and benign classes after review filtering
        if label == 1:
            high_review_pathogenic_variants += 1
        else:
            high_review_benign_variants += 1


print("CLINVAR BENCHMARK FILTERING SUMMARY")
print("-----------------------------------")
print("Total ClinVar variants:", total_variants)
print("Missense variants:", missense_variants)
print("Missense SNVs:", snv_variants)
print("Valid one-base REF/ALT SNVs:", allele_valid_variants)
print("Accepted clinical significance:", clinical_significance_variants)
print("2-4-star variants:", high_review_variants)

print()
print("CLASS COUNTS BEFORE REVIEW-STATUS FILTER")
print("----------------------------------------")
print("Pathogenic / likely pathogenic:", pathogenic_variants)
print("Benign / likely benign:", benign_variants)

print()
print("CLASS COUNTS AFTER 2-4-STAR REVIEW-STATUS FILTER")
print("------------------------------------------------")
print(
    "Pathogenic / likely pathogenic:",
    high_review_pathogenic_variants
)
print(
    "Benign / likely benign:",
    high_review_benign_variants
)

class_count_verification = (
    high_review_pathogenic_variants
    + high_review_benign_variants
)

print("Class-count verification:", class_count_verification)
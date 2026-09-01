import gzip


vcf_file = "data/raw/clinvar_20260822.vcf.gz"

total_variants = 0
missense_variants = 0
snv_variants = 0
one_base_variants = 0
canonical_allele_variants = 0
noncanonical_allele_variants = 0
clinical_significance_variants = 0
high_review_variants = 0

pathogenic_variants = 0
benign_variants = 0

strict_pathogenic_variants = 0
strict_benign_variants = 0

canonical_bases = {"A", "C", "G", "T"}

pathogenic_classes = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
}

benign_classes = {
    "Benign",
    "Likely_benign",
    "Benign/Likely_benign",
}

accepted_review_statuses = {
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
}


with gzip.open(vcf_file, "rt") as file:
    for line in file:

        # Skip metadata and column-header lines
        if line.startswith("#"):
            continue

        total_variants += 1

        # Split the VCF row into columns
        fields = line.rstrip("\n").split("\t")

        ref = fields[3]
        alt = fields[4]
        info = fields[7]

        # Parse the INFO field
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

        # Filter 3: REF and ALT must each contain exactly one character
        if len(ref) != 1 or len(alt) != 1:
            continue

        one_base_variants += 1

        # Filter 4: REF and ALT must be canonical DNA nucleotides
        #
        # This removes records containing N or ".". In particular, ALT="."
        # represents no alternate allele and cannot define a nucleotide
        # substitution for the missense-SNV benchmark.
        if ref not in canonical_bases or alt not in canonical_bases:
            noncanonical_allele_variants += 1
            continue

        canonical_allele_variants += 1

        # Filter 5: accepted binary clinical significance
        if clinical_significance in pathogenic_classes:
            label = 1
            pathogenic_variants += 1

        elif clinical_significance in benign_classes:
            label = 0
            benign_variants += 1

        else:
            continue

        clinical_significance_variants += 1

        # Filter 6: strict 2-4-star review status
        if review_status not in accepted_review_statuses:
            continue

        high_review_variants += 1

        if label == 1:
            strict_pathogenic_variants += 1
        else:
            strict_benign_variants += 1


print("CLINVAR BENCHMARK FILTERING SUMMARY")
print("-----------------------------------")
print("Total ClinVar variants:", total_variants)
print("Missense variants:", missense_variants)
print("Missense SNVs:", snv_variants)
print("One-character REF/ALT SNVs:", one_base_variants)
print("Excluded noncanonical one-character alleles:", noncanonical_allele_variants)
print("Canonical A/C/G/T missense SNVs:", canonical_allele_variants)
print("Accepted clinical significance:", clinical_significance_variants)
print("Strict 2-4-star variants:", high_review_variants)

print("\nCLASS COUNTS BEFORE REVIEW-STATUS FILTER")
print("----------------------------------------")
print("Pathogenic / likely pathogenic:", pathogenic_variants)
print("Benign / likely benign:", benign_variants)
print(
    "Class-count verification:",
    pathogenic_variants + benign_variants,
)

print("\nCLASS COUNTS AFTER 2-4-STAR REVIEW-STATUS FILTER")
print("------------------------------------------------")
print("Pathogenic / likely pathogenic:", strict_pathogenic_variants)
print("Benign / likely benign:", strict_benign_variants)
print(
    "Class-count verification:",
    strict_pathogenic_variants + strict_benign_variants,
)
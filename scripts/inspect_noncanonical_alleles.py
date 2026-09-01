import gzip
from collections import Counter


vcf_file = "data/raw/clinvar_20260822.vcf.gz"

expected_strict_total = 65225

valid_nucleotides = {"A", "C", "G", "T"}

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

accepted_review_statuses = {
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline"
}


all_noncanonical_records = []
strict_noncanonical_records = []

strict_candidate_count = 0


def parse_info_field(info):
    """Convert a VCF INFO string into a dictionary."""

    info_dict = {}

    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = True

    return info_dict


with gzip.open(vcf_file, "rt") as file:

    for line in file:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 8:
            continue

        chrom = fields[0]
        pos = fields[1]
        clinvar_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        info = fields[7]

        info_dict = parse_info_field(info)

        molecular_consequence = info_dict.get("MC", "")
        variant_type = info_dict.get("CLNVC", "")
        clinical_significance = info_dict.get("CLNSIG", "")
        review_status = info_dict.get("CLNREVSTAT", "")
        gene_info = info_dict.get("GENEINFO", "")

        # Reproduce Filters 1-3
        if "missense_variant" not in molecular_consequence:
            continue

        if variant_type != "single_nucleotide_variant":
            continue

        if len(ref) != 1 or len(alt) != 1:
            continue

        is_noncanonical = (
            ref not in valid_nucleotides
            or alt not in valid_nucleotides
        )

        record = {
            "chrom": chrom,
            "pos": pos,
            "clinvar_id": clinvar_id,
            "ref": ref,
            "alt": alt,
            "clinical_significance": clinical_significance,
            "review_status": review_status,
            "gene_info": gene_info,
            "molecular_consequence": molecular_consequence
        }

        # Count noncanonical alleles before Filters 4 and 5
        if is_noncanonical:
            all_noncanonical_records.append(record)

        # Filter 4: binary clinical significance
        if clinical_significance in pathogenic_classes:
            label = "pathogenic"

        elif clinical_significance in benign_classes:
            label = "benign"

        else:
            continue

        # Filter 5: strict review status
        if review_status not in accepted_review_statuses:
            continue

        strict_candidate_count += 1

        # Determine whether a strict candidate is noncanonical
        if is_noncanonical:
            record["label"] = label
            strict_noncanonical_records.append(record)


if strict_candidate_count != expected_strict_total:
    raise RuntimeError(
        "Strict benchmark count changed: "
        f"expected {expected_strict_total}, "
        f"observed {strict_candidate_count}"
    )


all_allele_patterns = Counter(
    (record["ref"], record["alt"])
    for record in all_noncanonical_records
)

all_clnsig_values = Counter(
    record["clinical_significance"]
    for record in all_noncanonical_records
)

all_review_values = Counter(
    record["review_status"]
    for record in all_noncanonical_records
)

strict_allele_patterns = Counter(
    (record["ref"], record["alt"])
    for record in strict_noncanonical_records
)


print("NONCANONICAL ALLELE INSPECTION")
print("------------------------------")
print("Reproduced strict candidate count:", strict_candidate_count)
print(
    "Noncanonical records after Filters 1-3:",
    len(all_noncanonical_records)
)
print(
    "Noncanonical records in strict 2-4-star benchmark:",
    len(strict_noncanonical_records)
)

print()
print("ALL NONCANONICAL REF/ALT PATTERNS")
print("---------------------------------")

for (ref, alt), count in sorted(all_allele_patterns.items()):
    print(f"{ref}>{alt}: {count}")

print()
print("CLNSIG VALUES AMONG ALL NONCANONICAL RECORDS")
print("---------------------------------------------")

for value, count in all_clnsig_values.most_common():
    print(f"{value or '[missing]'}: {count}")

print()
print("CLNREVSTAT VALUES AMONG ALL NONCANONICAL RECORDS")
print("-------------------------------------------------")

for value, count in all_review_values.most_common():
    print(f"{value or '[missing]'}: {count}")

print()
print("STRICT NONCANONICAL REF/ALT PATTERNS")
print("------------------------------------")

if strict_allele_patterns:

    for (ref, alt), count in sorted(strict_allele_patterns.items()):
        print(f"{ref}>{alt}: {count}")

else:
    print("No noncanonical alleles occur in the strict benchmark.")

print()
print("STRICT NONCANONICAL RECORD DETAILS")
print("----------------------------------")

if strict_noncanonical_records:

    for record in strict_noncanonical_records:

        print(
            "\t".join(
                [
                    record["chrom"],
                    record["pos"],
                    record["clinvar_id"],
                    record["ref"],
                    record["alt"],
                    record["label"],
                    record["clinical_significance"],
                    record["review_status"],
                    record["gene_info"],
                    record["molecular_consequence"]
                ]
            )
        )

else:
    print("No strict noncanonical records require exclusion.")

print()
print("INSPECTION COMPLETE")
print("-------------------")
print("PASS: The strict benchmark total was reproduced.")
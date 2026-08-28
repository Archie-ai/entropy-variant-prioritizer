import gzip

vcf_file = "data/raw/clinvar_20260822.vcf.gz"

with gzip.open(vcf_file, "rt") as file:

    for line in file:

        # Skip metadata lines
        if line.startswith("##"):
            continue

        # Find and print the VCF column header
        if line.startswith("#CHROM"):
            print("VCF COLUMNS:")
            print(line.strip())
            continue

        # Split the first variant row into its 8 VCF columns
        fields = line.strip().split("\t")

        chrom = fields[0]
        pos = fields[1]
        clinvar_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        qual = fields[5]
        filter_status = fields[6]
        info = fields[7]

        # Print each column separately
        print("\nFIRST CLINVAR VARIANT")
        print("---------------------")
        print("Chromosome:", chrom)
        print("Position:", pos)
        print("ClinVar ID:", clinvar_id)
        print("Reference allele:", ref)
        print("Alternate allele:", alt)
        print("Quality:", qual)
        print("Filter:", filter_status)
        print("INFO:", info)

        # Stop after the first variant
        break
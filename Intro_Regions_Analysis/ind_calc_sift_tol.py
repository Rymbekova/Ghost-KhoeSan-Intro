import subprocess
import os
import pandas as pd
import numpy as np
import re
import sys

# Check if an individual name is provided as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <INDIVIDUAL>")
    sys.exit(1)

INDIVIDUAL = sys.argv[1]  # Get the individual name from command-line arguments

# Update BED file paths dynamically based on the individual
REAL_BED_FILE = f"/lisc/scratch/admixlab/aigerim/sstar-analysis-main/African/San21/{INDIVIDUAL}/{INDIVIDUAL}.all.merged.tracts.bed"
RANDOM_BED_FILE = f"/lisc/scratch/admixlab/aigerim/sstar-analysis-main/African/San21/{INDIVIDUAL}/{INDIVIDUAL}.random.merged.tracts.bed"
VCF_DIR = "/lisc/scratch/admixlab/aigerim/ensembl-vep/vcf"
VCF_FILE_TEMPLATE = "sstar_subset_filter1_{chr}.vcf.gz"
CHROMOSOMES = [str(i) for i in range(1, 23)]
SIFT_BED_FILE = "SIFT_tol_sites.txt"
GENOME_FILE = "genome.txt"  # Update with the path to genome.txt

def extract_genotypes(bed_file):
    genotype_data = []
    for chr in CHROMOSOMES:
        print(f"Processing chromosome {chr}")
        vcf_file = os.path.join(VCF_DIR, VCF_FILE_TEMPLATE.format(chr=chr))

        chr_real_bed = f"chr{chr}_regions.bed"
        chr_sift_bed = f"chr{chr}_sift.bed"
        chr_intersect_bed_with_allele = f"chr{chr}_intersect_with_allele.bed"

        # Extract chromosome-specific regions
        with open(bed_file, 'r') as bed_in, open(chr_real_bed, 'w') as bed_out:
            for line in bed_in:
                if line.startswith(f"chr{chr}\t") or line.startswith(f"{chr}\t"):
                    bed_out.write(line)
        with open(SIFT_BED_FILE, 'r') as sift_in, open(chr_sift_bed, 'w') as sift_out:
            for line in sift_in:
                if line.startswith(f"chr{chr}\t") or line.startswith(f"{chr}\t"):
                    sift_out.write(line)

        # Sort the BED files
        sorted_chr_real_bed = f"sorted_chr{chr}_regions.bed"
        sorted_chr_sift_bed = f"sorted_chr{chr}_sift.bed"

        subprocess.run([
            "bedtools", "sort",
            "-i", chr_real_bed,
            "-g", GENOME_FILE
        ], stdout=open(sorted_chr_real_bed, 'w'))

        subprocess.run([
            "bedtools", "sort",
            "-i", chr_sift_bed,
            "-g", GENOME_FILE
        ], stdout=open(sorted_chr_sift_bed, 'w'))

        # Intersect the sorted BED files and include both entries
        chr_intersect_bed_with_allele = f"chr{chr}_intersect_with_allele.bed"
        subprocess.run([
            "bedtools", "intersect",
            "-a", sorted_chr_real_bed,
            "-b", sorted_chr_sift_bed,
            "-wa",
            "-wb",
            "-sorted"
        ], stdout=open(chr_intersect_bed_with_allele, 'w'))

        # Check if the intersected BED file has any content
        if os.path.getsize(chr_intersect_bed_with_allele) == 0:
            print(f"No overlapping regions found for chromosome {chr}. Skipping...")
            continue

        # Read the intersected BED file
        try:
            intersect_df = pd.read_csv(chr_intersect_bed_with_allele, sep='\t', header=None)
        except pd.errors.EmptyDataError:
            print(f"No data in intersected BED file for chromosome {chr}. Skipping...")
            continue

        intersect_df = intersect_df.drop_duplicates()

        # Save positions to a temporary file
        temp_bed = f"temp_chr{chr}_positions.txt"
        intersect_df[[3, 4, 5]].to_csv(temp_bed, sep='\t', index=False, header=False)

        # Extract genotypes using bcftools (include square brackets around %GT)
        temp_output = f"chr{chr}_genotypes.txt"
        cmd = [
            "bcftools", "view", "-m2", "-M2", "-v", "snps",
            "-R", temp_bed,
            "-s", INDIVIDUAL,
            vcf_file
        ]

        query_cmd = [
            "bcftools", "query",
            "-f", '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n'  # Define the output format
        ]

        with open(temp_output, 'w') as f_out:
            # Start the bcftools view process
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # Start the bcftools query process, taking input from p1's stdout
            p2 = subprocess.Popen(query_cmd, stdin=p1.stdout, stdout=f_out, stderr=subprocess.PIPE)
            # Close p1's stdout to allow p1 to receive a SIGPIPE if p2 exits
            p1.stdout.close()
            # Wait for both processes to complete
            _, err1 = p1.communicate()
            _, err2 = p2.communicate()


        # Check if the temp_output file was created successfully and is not empty
        if not os.path.exists(temp_output) or os.path.getsize(temp_output) == 0:
            print(f"No genotype data found for chromosome {chr}. Skipping...")
            continue

        # Read genotype data
        try:
            genotype_df = pd.read_csv(
                temp_output, sep='\t', names=['CHROM', 'POS', 'REF', 'ALT', 'GT']
            )
        except pd.errors.EmptyDataError:
            print(f"No genotype data available in file {temp_output}. Skipping...")
            continue

        # Ensure 'POS' is numeric and drop rows with NaN POS
        genotype_df['POS'] = pd.to_numeric(genotype_df['POS'], errors='coerce')
        genotype_df = genotype_df.dropna(subset=['POS'])
        genotype_df['POS'] = genotype_df['POS'].astype(int)
        #genotype_df['CHROM'] = genotype_df['CHROM'].str.replace('chr', '', regex=False)
        genotype_df = genotype_df.drop_duplicates(subset=['CHROM', 'POS'])

        # Debugging statements
        print(f"Chromosome {chr}: Number of intersected positions: {len(intersect_df)}")
        print(f"Chromosome {chr}: Genotype data count: {len(genotype_df)}")
        #if not genotype_df.empty:
            #print(f"Chromosome {chr}: Sample genotype data:")
            #print(genotype_df.head())

        genotype_data.append(genotype_df)

        # Clean up temporary files
        os.remove(chr_real_bed)
        os.remove(chr_sift_bed)
        os.remove(sorted_chr_real_bed)
        os.remove(sorted_chr_sift_bed)
        os.remove(chr_intersect_bed_with_allele)
        os.remove(temp_bed)
        os.remove(temp_output)

    return genotype_data

def process_genotypes(genotype_data, output_file):
    df = pd.concat(genotype_data, ignore_index=True)

    if df.empty:
        print("No genotype data found.")
        return 0

    # Ensure 'POS' is of integer type
    df['POS'] = df['POS'].astype(int)
    def count_deleterious_alleles(gt):
        if gt in ['./.', '.|.', '0/0', '0|0']:
            return 0
        return 1

    # Apply the function to the 'GT' column
    df['Recoded_GT'] = df['GT'].apply(count_deleterious_alleles)

    # Debugging statement
    print(f"Total number of genotypes processed: {len(df)}")
    if not df.empty:
        print(f"Sample recoded genotypes:")
        print(df[['CHROM', 'POS', 'GT', 'Recoded_GT']].head())

    # Save the DataFrame to a CSV file
    df.to_csv(output_file, index=False)

    # Calculate the sum of the recoded genotypes
    total_sum = df['Recoded_GT'].sum(skipna=True)
    return total_sum

# Extract and process genotypes for real regions
real_genotype_data = extract_genotypes(REAL_BED_FILE)
real_output_file = "sifl_tol_real_regions_genotypes.csv"
real_sum = process_genotypes(real_genotype_data, real_output_file)

# Extract and process genotypes for random regions
random_genotype_data = extract_genotypes(RANDOM_BED_FILE)
random_output_file = "sift_tol_random_regions_genotypes.csv"
random_sum = process_genotypes(random_genotype_data, random_output_file)

# Compare sums
print(f"Sum over real regions: {real_sum}")
print(f"Sum over random regions: {random_sum}")
print(f"Difference between real and random sums: {real_sum - random_sum}")

import subprocess
import os
import pandas as pd
import sys
import tempfile
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

# Check if an individual name is provided as a command-line argument
if len(sys.argv) != 2:
    logging.error("Usage: python script.py <INDIVIDUAL>")
    sys.exit(1)

INDIVIDUAL = sys.argv[1]  # Get the individual name from command-line arguments

# Update BED file paths dynamically based on the individual
REAL_BED_FILE = f"/lisc/scratch/admixlab/aigerim/sstar-analysis-main/African/San21/{INDIVIDUAL}/{INDIVIDUAL}.all.merged.tracts.bed"
RANDOM_BED_FILE = f"/lisc/scratch/admixlab/aigerim/sstar-analysis-main/African/San21/{INDIVIDUAL}/{INDIVIDUAL}.random.merged.tracts.bed"
VCF_DIR = "/lisc/scratch/admixlab/aigerim/ensembl-vep/vcf"
VCF_FILE_TEMPLATE = "sstar_subset_filter1_{chr}.vcf.gz"
CHROMOSOMES = [str(i) for i in range(1, 23)]
HIGH_BED_FILE = "../African/final_gerp_high_auto.txt"
GENOME_FILE = "genome.txt"  # Update with the path to genome.txt
NUM_WORKERS = 8  # Set the number of parallel workers to 8

def count_deleterious_alleles(gt):
    if gt in ['./.', '.|.', '0/0', '0|0']:
        return 0
    return 1

def process_chromosome(chr, bed_file, individual, high_bed_file, vcf_dir, vcf_file_template, genome_file, temp_dir):
    """
    Processes a single chromosome: extracts regions, intersects BED files, extracts genotypes, and recodes them.
    Returns a pandas DataFrame with the genotype data or None if no data is found.
    """
    CHROMOSOME = str(chr)
    VCF_FILE = os.path.join(vcf_dir, vcf_file_template.format(chr=CHROMOSOME))
    
    # Define temporary file names within the temporary directory
    chr_real_bed = os.path.join(temp_dir, f"chr{CHROMOSOME}_regions.bed")
    chr_high_bed = os.path.join(temp_dir, f"chr{CHROMOSOME}_high.bed")
    sorted_chr_real_bed = os.path.join(temp_dir, f"sorted_chr{CHROMOSOME}_regions.bed")
    sorted_chr_high_bed = os.path.join(temp_dir, f"sorted_chr{CHROMOSOME}_high.bed")
    chr_intersect_bed_with_allele = os.path.join(temp_dir, f"chr{CHROMOSOME}_intersect_with_allele.bed")
    temp_bed = os.path.join(temp_dir, f"temp_chr{CHROMOSOME}_positions.txt")
    temp_output = os.path.join(temp_dir, f"chr{CHROMOSOME}_genotypes.txt")

    # Extract chromosome-specific regions from the real BED file
    with open(bed_file, 'r') as bed_in, open(chr_real_bed, 'w') as bed_out:
        for line in bed_in:
            if line.startswith(f"chr{CHROMOSOME}\t") or line.startswith(f"{CHROMOSOME}\t"):
                bed_out.write(line)
    
    # Extract chromosome-specific regions from the HIGH BED file
    with open(high_bed_file, 'r') as high_in, open(chr_high_bed, 'w') as high_out:
        for line in high_in:
            if line.startswith(f"chr{CHROMOSOME}\t") or line.startswith(f"{CHROMOSOME}\t"):
                high_out.write(line)

    # Sort the BED files using bedtools
    subprocess.run(
        ["bedtools", "sort", "-i", chr_real_bed, "-g", genome_file],
        stdout=open(sorted_chr_real_bed, 'w'),
        check=True
    )

    subprocess.run(
        ["bedtools", "sort", "-i", chr_high_bed, "-g", genome_file],
        stdout=open(sorted_chr_high_bed, 'w'),
        check=True
    )

    # Intersect the sorted BED files
    subprocess.run(
        [
            "bedtools", "intersect",
            "-a", sorted_chr_real_bed,
            "-b", sorted_chr_high_bed,
            "-wa", "-wb", "-sorted"
        ],
        stdout=open(chr_intersect_bed_with_allele, 'w'),
        check=True
    )

    # Check if the intersected BED file has any content
    if os.path.getsize(chr_intersect_bed_with_allele) == 0:
        logging.info(f"No overlapping regions found for chromosome {CHROMOSOME}. Skipping...")
        return None

    # Read the intersected BED file
    try:
        intersect_df = pd.read_csv(chr_intersect_bed_with_allele, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        logging.info(f"No data in intersected BED file for chromosome {CHROMOSOME}. Skipping...")
        return None

    intersect_df = intersect_df.drop_duplicates()

    # Save positions to a temporary file
    intersect_df[[3, 4, 5]].to_csv(temp_bed, sep='\t', index=False, header=False)

    # Extract genotypes using bcftools
    cmd = [
        "bcftools", "view", "-m2", "-M2", "-v", "snps",
        "-R", temp_bed,
        "-s", individual,
        VCF_FILE
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
        logging.info(f"No genotype data found for chromosome {CHROMOSOME}. Skipping...")
        return None

    # Read genotype data
    try:
        genotype_df = pd.read_csv(
            temp_output, sep='\t', names=['CHROM', 'POS', 'REF', 'ALT', 'GT']
        )
    except pd.errors.EmptyDataError:
        logging.info(f"No genotype data available in file {temp_output}. Skipping...")
        return None

    # Ensure 'POS' is numeric and drop rows with NaN POS
    genotype_df['POS'] = pd.to_numeric(genotype_df['POS'], errors='coerce')
    genotype_df = genotype_df.dropna(subset=['POS'])
    genotype_df['POS'] = genotype_df['POS'].astype(int)
    genotype_df = genotype_df.drop_duplicates(subset=['CHROM', 'POS'])

    # Recoding genotypes
    genotype_df['Recoded_GT'] = genotype_df['GT'].apply(count_deleterious_alleles)

    return genotype_df

def extract_genotypes_parallel(bed_file, individual, high_bed_file, vcf_dir, vcf_file_template, genome_file, max_workers=8):
    """
    Extracts genotypes in parallel across chromosomes using a specified number of worker processes.
    Returns a list of pandas DataFrames containing genotype data.
    """
    genotype_data = []

    with tempfile.TemporaryDirectory() as temp_dir:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all chromosome tasks
            futures = {
                executor.submit(
                    process_chromosome,
                    chr,
                    bed_file,
                    individual,
                    high_bed_file,
                    vcf_dir,
                    vcf_file_template,
                    genome_file,
                    temp_dir
                ): chr for chr in CHROMOSOMES
            }

            for future in as_completed(futures):
                chr = futures[future]
                try:
                    result = future.result()
                    if result is not None:
                        genotype_data.append(result)
                    else:
                        logging.info(f"Chromosome {chr} returned no data.")
                except Exception as e:
                    logging.error(f"Error processing chromosome {chr}: {e}")

    return genotype_data

def process_genotypes(genotype_data, output_file):
    """
    Processes the collected genotype data: concatenates DataFrames, recodes genotypes, and calculates the sum.
    Saves the processed data to a CSV file.
    Returns the total sum of recoded genotypes.
    """
    if not genotype_data:
        logging.warning("No genotype data found.")
        return 0

    df = pd.concat(genotype_data, ignore_index=True)

    if df.empty:
        logging.warning("No genotype data found after concatenation.")
        return 0

    # Ensure 'POS' is of integer type
    df['POS'] = df['POS'].astype(int)

    # Recoding genotypes
    df['Recoded_GT'] = df['GT'].apply(count_deleterious_alleles)

    # Debugging statements
    logging.info(f"Total number of genotypes processed: {len(df)}")
    if not df.empty:
        logging.info("Sample recoded genotypes:")
        logging.info(df[['CHROM', 'POS', 'GT', 'Recoded_GT']].head())

    # Save the DataFrame to a CSV file with the individual's name as a prefix
    df.to_csv(output_file, index=False)

    # Calculate the sum of the recoded genotypes
    total_sum = df['Recoded_GT'].sum(skipna=True)
    return total_sum

def main():
    # Set the number of parallel workers to 8
    num_workers = NUM_WORKERS
    logging.info(f"Setting number of workers to: {num_workers}")

    # Extract and process genotypes for real regions in parallel
    logging.info("Processing real regions in parallel...")
    real_genotype_data = extract_genotypes_parallel(
        REAL_BED_FILE,
        INDIVIDUAL,
        HIGH_BED_FILE,
        VCF_DIR,
        VCF_FILE_TEMPLATE,
        GENOME_FILE,
        max_workers=num_workers
    )
    real_output_file = f"{INDIVIDUAL}_gerp_high_real_regions_genotypes.csv"
    real_sum = process_genotypes(real_genotype_data, real_output_file)

    # Extract and process genotypes for random regions in parallel
    logging.info("Processing random regions in parallel...")
    random_genotype_data = extract_genotypes_parallel(
        RANDOM_BED_FILE,
        INDIVIDUAL,
        HIGH_BED_FILE,
        VCF_DIR,
        VCF_FILE_TEMPLATE,
        GENOME_FILE,
        max_workers=num_workers
    )
    random_output_file = f"{INDIVIDUAL}_gerp_high_random_regions_genotypes.csv"
    random_sum = process_genotypes(random_genotype_data, random_output_file)

    # Compare sums
    logging.info(f"Sum over real regions: {real_sum}")
    logging.info(f"Sum over random regions: {random_sum}")
    logging.info(f"Difference between real and random sums: {real_sum - random_sum}")

if __name__ == "__main__":
    main()

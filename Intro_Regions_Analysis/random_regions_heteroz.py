import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import os
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import subprocess
import csv
import argparse

# Configure logging to show detailed messages
logging.basicConfig(
    level=logging.INFO,  # Set to INFO or WARNING to reduce verbosity
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)

def load_chrom_limits(chrom_limits_file):
    """
    Load chromosomal limits from a file into a dictionary.
    """
    try:
        limits_df = pd.read_csv(
            chrom_limits_file,
            sep='\t',
            header=None,
            names=['chrom', 'first_pos', 'last_pos']
        )
        chrom_limits = limits_df.set_index('chrom').T.to_dict('list')
        logging.info(f"Chromosomal limits loaded from {chrom_limits_file}")
    except Exception as e:
        logging.error(f"Error loading chromosomal limits: {e}")
        raise e
    return chrom_limits

def load_informative_snps(snp_dir, chromosomes):
    """
    Load informative SNP positions for each chromosome from pre-parsed SNP files into a dictionary.
    Each SNP file should be named like 'chr1_informative_snps.txt', containing one SNP position per line.
    """
    snp_dict = {}
    for chrom in chromosomes:
        chrom_num = chrom.replace('chr', '')
        snp_file = os.path.join(snp_dir, f"chr{chrom_num}_informative_snps.txt")
        
        if not os.path.exists(snp_file):
            logging.warning(f"SNP file {snp_file} does not exist. Skipping chromosome {chrom}.")
            snp_dict[chrom] = np.array([])
            continue

        try:
            snps = np.loadtxt(snp_file, usecols=1, dtype=int)
            snp_dict[chrom] = np.array(sorted(snps))
            logging.info(f"Loaded {len(snps)} informative SNPs for {chrom} from {snp_file}")
        except Exception as e:
            logging.error(f"Error loading SNPs from {snp_file}: {e}")
            snp_dict[chrom] = np.array([])
    return snp_dict

def count_informative_snps_cached(snp_positions, start, end):
    """
    Count the number of informative SNPs within a specified region using preloaded SNP positions.
    """
    if len(snp_positions) == 0:
        return 0
    left = np.searchsorted(snp_positions, start, side='left')
    right = np.searchsorted(snp_positions, end, side='right')
    return right - left

def calculate_heterozygosity_real(region, vcf_file, individual):
    """
    Calculate heterozygosity for a real region and individual using heterozygosity_real.py.
    """
    command = ["python", "heterozygosity_real.py", region, vcf_file, individual]
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode == 0:
        output = result.stdout.strip().split()
        if output:
            return float(output[-1])
    return "Error in calculate_heterozygosity_real"


def calculate_heterozygosity_random(region, vcf_file, individual):
    """
    Calculate heterozygosity for a random region and individual using heterozygosity_real.py.
    """
    command = ["python", "heterozygosity_random.py", region, vcf_file, individual]
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode == 0:
        output = result.stdout.strip().split()
        if output:
            return float(output[-1])
    return "Error in calculate_heterozygosity_random"


def generate_random_start(row, snp_positions, chrom_interval_tree, real_snp_count, real_regions, snp_tol=0.05, max_attempts_per_iteration=10000, num_iterations=100):
    """
    Generate random regions for a genomic region over a specified number of iterations (default 100).
    Ensures no overlap with existing intervals (including real regions) and the SNP count is within ±5% of the real region's SNP count.
    If max attempts are reached, reset the interval tree and keep retrying until a valid random region is found.
    """
    range_length = row['original_length']
    min_start, max_start = row['first_pos'], row['last_pos'] - range_length
    lower_bound = int(real_snp_count * (1 - snp_tol))
    upper_bound = int(real_snp_count * (1 + snp_tol))
    
    generated_regions = []  # Store generated random regions and their SNP counts

    # Log starting information
    logging.info(f"Starting random region generation for {row['chrom']} with {num_iterations} iterations and {max_attempts_per_iteration} attempts per iteration.")
    
    for i in range(num_iterations):  # Iterate num_iterations times to generate random regions
        while True:  # Continue trying until a valid random region is found
            attempts = 0
            
            while attempts < max_attempts_per_iteration:
                random_start = np.random.randint(min_start, max_start + 1)
                random_end = random_start + range_length

                # Check for overlap with existing intervals
                if chrom_interval_tree.overlap(random_start, random_end):
                    attempts += 1
                    continue

                # Count SNPs in the random region
                snp_count = count_informative_snps_cached(snp_positions, random_start, random_end)

                # Check if SNP count is within the tolerance
                if lower_bound <= snp_count <= upper_bound:
                    chrom_interval_tree.add(Interval(random_start, random_end))  # Add interval to the tree to prevent overlaps
                    generated_regions.append((random_start, snp_count))  # Store the random region and SNP count
                    logging.info(f"Random region {i+1}: {random_start}-{random_end} (SNP count: {snp_count})")
                    break  # Exit the while loop and move to the next iteration if a valid region is found

                attempts += 1
            
            # If max attempts are reached, reset the interval tree and retry until solution found
            if attempts == max_attempts_per_iteration:
                logging.warning(f"No valid random region found for {row['chrom']} after {max_attempts_per_iteration} attempts in iteration {i+1}. Resetting interval tree and retrying.")

                # Reset the interval tree (clear random intervals)
                chrom_interval_tree.clear()

                # Re-add the real regions to the tree to prevent overlap with them
                for real_start, real_end in real_regions:
                    chrom_interval_tree.add(Interval(real_start, real_end))
                    logging.info(f"Re-added real region: {real_start}-{real_end}.")
            
            # If a valid region was found, exit the infinite while loop and continue to the next iteration
            else:
                break

    return generated_regions

def run_script(script_name, region, vcf_dir, chrom, individual):
    """
    Runs the given Python script and returns the heterozygosity value or an error.
    """
    vcf_file = os.path.join(vcf_dir, f"25KS.48RHG.74comp.HCBP.{chrom.replace('chr', '')}.recalSNP99.9.recalINDEL99.0.vcf.gz")
    command = ["python", script_name, region, vcf_file, individual]
    
    # Run the subprocess and capture the output
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode == 0 and result.stdout.strip():
        output = result.stdout.strip()
        if output.replace('.', '', 1).isdigit():  # Check if the output is a valid float
            return float(output)
        logging.error(f"Failed to parse heterozygosity value from {script_name} output: {output}")
    else:
        logging.error(f"Error running {script_name} for region {region}: {result.stderr}")

    return f"Error in run_script for {script_name}"

def save_to_csv(output_file, output_rows, header):
    """
    Saves the given data (output_rows) to a CSV file with the specified header.
    """
    try:
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(header)  # Write the header
            writer.writerows(output_rows)  # Write the data rows
        print(f"Data successfully saved to {output_file}")
    except Exception as e:
        print(f"Error saving to CSV file: {e}")

def process_individual(individual, snp_dict, snp_tol, snp_dir, chrom_limits, chromosomes):
    input_file = os.path.join(individual, f"{individual}.all.merged.tracts.bed")
    output_file = os.path.join(individual, f"{individual}_heterozygosity_reset.csv")
    vcf_dir = '/lisc/scratch/admixlab/aigerim/African/'

    try:
        # Load the regions for the individual
        ranges_df = pd.read_csv(input_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
        ranges_df['original_length'] = ranges_df['end'] - ranges_df['start']
        ranges_df[['first_pos', 'last_pos']] = ranges_df['chrom'].map(
            lambda x: chrom_limits.get(x, [np.nan, np.nan])
        ).tolist()
    except Exception as e:
        logging.error(f"Error loading/processing tracts for {individual}: {e}")
        return

    # Prepare interval trees for each chromosome to track regions
    chrom_interval_trees = {chrom: IntervalTree() for chrom in chromosomes}
    
    # Add real regions to the interval tree
    for chrom, start, end in zip(ranges_df['chrom'], ranges_df['start'], ranges_df['end']):
        chrom_interval_trees[chrom].addi(start, end)

    output_rows = []

    for chrom in chromosomes:
        snp_positions = snp_dict.get(chrom, np.array([]))
        if snp_positions.size == 0 or ranges_df[ranges_df['chrom'] == chrom].empty:
            continue

        # Filter ranges for the current chromosome
        chrom_ranges = ranges_df[ranges_df['chrom'] == chrom].copy()
        chrom_ranges['snp_count_real'] = chrom_ranges.apply(
            lambda row: count_informative_snps_cached(snp_positions, row['start'], row['end']),
            axis=1
        )

        # Extract real regions for this chromosome to pass to generate_random_start
        real_regions = [(row['start'], row['end']) for _, row in chrom_ranges.iterrows()]

        for idx, row in chrom_ranges.iterrows():
            real_region = f"{row['chrom']}:{row['start']}-{row['end']}"
            snp_count_real = row['snp_count_real']
            real_het = run_script("heterozygosity_real.py", real_region, vcf_dir, chrom, individual)
            logging.info(f"Real heterozygosity for {real_region} in {individual}: {real_het}")

            # Store real heterozygosity and random heterozygosity values
            random_data = generate_random_start(row, snp_positions, chrom_interval_trees[chrom], snp_count_real, real_regions, snp_tol, 100)
            random_het_values = []
            
            for start, _ in random_data:
                if not np.isnan(start):
                    random_region = f"{row['chrom']}:{int(start)}-{int(start + row['original_length'])}"
                    random_het = run_script("heterozygosity_random.py", random_region, vcf_dir, chrom, individual)
                    logging.info(f"Random heterozygosity for {random_region} in {individual}: {random_het}")
                    random_het_values.append(random_het)

            # Append real and random heterozygosity to output_rows
            output_rows.append([real_region, real_het] + random_het_values)

    # Save all results for this individual to a CSV file
    save_to_csv(output_file, output_rows, ['Real Region', 'Real Heterozygosity'] + [f'Random_Het_{i+1}' for i in range(100)])


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process an individual for heterozygosity calculation.")
    parser.add_argument('individual', type=str, help="The name of the individual to process (e.g., KSP092).")

    # Parse the arguments
    args = parser.parse_args()

    chrom_limits_file = 'chrom_limits.txt'
    snp_dir = 'SNP_files'
    chromosomes = [f'chr{i}' for i in range(1, 23)]
    snp_tolerance = 0.05

    try:
        chrom_limits = load_chrom_limits(chrom_limits_file)
    except Exception:
        logging.error("Failed to load chromosomal limits. Exiting.")
        return

    snp_dict = load_informative_snps(snp_dir, chromosomes)

    # Call process_individual with the parsed individual name
    process_individual(
        individual=args.individual,
        snp_dict=snp_dict,
        snp_tol=snp_tolerance,
        snp_dir=snp_dir,
        chrom_limits=chrom_limits,
        chromosomes=chromosomes
    )

    logging.info(f"Completed processing for {args.individual}")

if __name__ == "__main__":
    main()

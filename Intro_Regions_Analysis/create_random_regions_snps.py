import pandas as pd
import numpy as np
import os
import argparse
import logging
from intervaltree import Interval, IntervalTree

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler()]
)

def load_chrom_limits(chrom_limits_file):
    """Load chromosomal limits from a file into a dictionary."""
    try:
        limits_df = pd.read_csv(
            chrom_limits_file, sep='\t', header=None, names=['chrom', 'first_pos', 'last_pos']
        )
        chrom_limits = limits_df.set_index('chrom').T.to_dict('list')
        logging.info(f"Chromosomal limits loaded from {chrom_limits_file}")
    except Exception as e:
        logging.error(f"Error loading chromosomal limits: {e}")
        raise e
    return chrom_limits

def load_informative_snps(snp_dir, chromosomes):
    """Load SNP positions for each chromosome from pre-parsed SNP files."""
    snp_dict = {}
    for chrom in chromosomes:
        snp_file = os.path.join(snp_dir, f"{chrom}_informative_snps.txt")
        
        if not os.path.exists(snp_file):
            logging.warning(f"SNP file {snp_file} does not exist. Skipping chromosome {chrom}.")
            snp_dict[chrom] = np.array([])
            continue

        try:
            snps = np.loadtxt(snp_file, usecols=1, dtype=int)
            snp_dict[chrom] = np.array(sorted(snps))
            logging.info(f"Loaded {len(snps)} SNPs for {chrom} from {snp_file}")
        except Exception as e:
            logging.error(f"Error loading SNPs from {snp_file}: {e}")
            snp_dict[chrom] = np.array([])
    return snp_dict

def count_informative_snps(snp_positions, start, end):
    """Count the number of SNPs within a specified region."""
    if len(snp_positions) == 0:
        return 0
    left = np.searchsorted(snp_positions, start, side='left')
    right = np.searchsorted(snp_positions, end, side='right')
    return right - left

def generate_random_region(chrom, region_length, target_snp_count, snp_positions, chrom_limits, interval_tree, real_regions, snp_tol=0.05, max_attempts=10000, max_retries=3):
    """
    Generate one random region with the same length and similar SNP count as the real region.
    If unable to find a suitable random region within max_attempts, reset the interval tree and retry.
    """
    min_start, max_start = chrom_limits[chrom][0], chrom_limits[chrom][1] - region_length
    lower_bound = int(target_snp_count * (1 - snp_tol))
    upper_bound = int(target_snp_count * (1 + snp_tol))
    
    for retry in range(max_retries):
        attempts = 0
        while attempts < max_attempts:
            random_start = np.random.randint(min_start, max_start + 1)
            random_end = random_start + region_length

            # Check for overlap with existing intervals
            if interval_tree.overlap(random_start, random_end):
                attempts += 1
                continue

            # Count SNPs in the random region
            snp_count = count_informative_snps(snp_positions, random_start, random_end)

            # Check if SNP count is within the tolerance
            if lower_bound <= snp_count <= upper_bound:
                interval_tree.add(Interval(random_start, random_end))
                return random_start, random_end

            attempts += 1

        # If max attempts reached, reset the interval tree and retry
        logging.warning(f"Unable to find a suitable random region for {chrom} within {max_attempts} attempts. Resetting interval tree and retrying (retry {retry + 1}/{max_retries}).")
        interval_tree.clear()
        
        # Re-add the real regions to the interval tree to prevent overlap
        for real_start, real_end in real_regions:
            interval_tree.add(Interval(real_start, real_end))
    
    # If all retries failed, raise an error
    raise ValueError(f"Unable to find a suitable random region for {chrom} after {max_retries} retries.")

def load_individuals(individuals_file):
    """Load individual names from a file."""
    with open(individuals_file, 'r') as file:
        individuals = [line.strip() for line in file if line.strip()]
    logging.info(f"Loaded {len(individuals)} individuals from {individuals_file}")
    return individuals

def main(individuals_file):
    chrom_limits_file = "chrom_limits.txt"
    snp_dir = "SNP_files"
    snp_tol = 0

    # Load chromosomal limits and SNP data
    chrom_limits = load_chrom_limits(chrom_limits_file)
    chromosomes = [f'chr{i}' for i in range(1, 23)]
    snp_dict = load_informative_snps(snp_dir, chromosomes)

    # Load individual names from file
    individuals = load_individuals(individuals_file)

    for individual in individuals:
        # Define paths for input and output files for each individual
        real_regions_file = f"{individual}/{individual}.all.merged.tracts.bed"
        output_bed_file = f"{individual}/{individual}.random.merged.tracts.bed"
        
        # Load real regions for the individual
        try:
            real_regions_df = pd.read_csv(real_regions_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
        except FileNotFoundError:
            logging.warning(f"Real regions file not found for {individual}. Skipping.")
            continue

        real_regions_df['length'] = real_regions_df['end'] - real_regions_df['start']
        bed_rows = []

        for chrom in chromosomes:
            chrom_regions = real_regions_df[real_regions_df['chrom'] == chrom]
            
            # Collect real regions for the interval tree
            real_regions = [(row['start'], row['end']) for _, row in chrom_regions.iterrows()]
            
            # Create an interval tree to track existing regions for overlap checking
            interval_tree = IntervalTree(Interval(start, end) for start, end in real_regions)

            for idx, row in chrom_regions.iterrows():
                start = row['start']
                end = row['end']
                region_length = row['length']
                
                if chrom not in snp_dict:
                    logging.warning(f"SNP data for {chrom} not found. Skipping region {chrom}:{start}-{end}")
                    continue

                # Count SNPs in the real region
                target_snp_count = count_informative_snps(snp_dict[chrom], start, end)

                # Generate a random region with similar SNP count and length
                try:
                    random_start, random_end = generate_random_region(
                        chrom, region_length, target_snp_count, snp_dict[chrom], chrom_limits, 
                        interval_tree=interval_tree, real_regions=real_regions, snp_tol=snp_tol
                    )
                    
                    if random_start is not None and random_end is not None:
                        bed_rows.append([chrom, random_start, random_end])
                    else:
                        logging.warning(f"Skipping {chrom}:{start}-{end} as no suitable random region was found.")

                except ValueError as e:
                    logging.warning(f"Failed to generate random region for {chrom}:{start}-{end}: {e}")

        # Save output in BED format for the individual
        bed_df = pd.DataFrame(bed_rows, columns=['chrom', 'start', 'end'])
        bed_df.to_csv(output_bed_file, sep='\t', header=False, index=False)
        logging.info(f"Random regions in BED format saved to {output_bed_file} for {individual}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate one set of random regions with similar length and SNP count to real regions for multiple individuals.")
    parser.add_argument("individuals_file", type=str, help="Path to the file containing individual names (one per line).")
    
    args = parser.parse_args()
    main(args.individuals_file)

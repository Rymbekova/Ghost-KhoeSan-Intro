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

def generate_random_region(chrom, region_length, target_snp_count, snp_positions,
                           chrom_limits, interval_tree, real_regions,
                           snp_tol=0.05, max_attempts=10000, max_retries=3):
    """
    Generate one random region with same length and similar SNP count.
    """
    min_start, max_start = chrom_limits[chrom][0], chrom_limits[chrom][1] - region_length
    lower_bound = int(target_snp_count * (1 - snp_tol))
    upper_bound = int(target_snp_count * (1 + snp_tol))

    for retry in range(max_retries):
        attempts = 0
        while attempts < max_attempts:
            random_start = np.random.randint(min_start, max_start + 1)
            random_end = random_start + region_length
            if interval_tree.overlap(random_start, random_end):
                attempts += 1
                continue
            snp_count = count_informative_snps(snp_positions, random_start, random_end)
            if lower_bound <= snp_count <= upper_bound:
                interval_tree.add(Interval(random_start, random_end))
                return random_start, random_end
            attempts += 1

        logging.warning(
            f"Failed to find region for {chrom} in {max_attempts} attempts. "
            f"Resetting interval tree and retrying (retry {retry+1}/{max_retries})."
        )
        interval_tree.clear()
        for start, end in real_regions:
            interval_tree.add(Interval(start, end))

    raise ValueError(f"Unable to find suitable random region for {chrom} after {max_retries} retries.")


def load_individuals(individuals_file):
    """Load individual names from a file."""
    with open(individuals_file, 'r') as f:
        individuals = [line.strip() for line in f if line.strip()]
    logging.info(f"Loaded {len(individuals)} individuals from {individuals_file}")
    return individuals


def main(individuals_file, num_replicates=100):
    chrom_limits_file = "chrom_limits.txt"
    snp_dir = "SNP_files"
    snp_tol = 0

    # Load limits and SNP data
    chrom_limits = load_chrom_limits(chrom_limits_file)
    chromosomes = [f'chr{i}' for i in range(1, 23)]
    snp_dict = load_informative_snps(snp_dir, chromosomes)
    individuals = load_individuals(individuals_file)

    for individual in individuals:
        real_file = f"{individual}/{individual}_overlap.bed"
        try:
            real_df = pd.read_csv(real_file, sep='\t', header=None,
                                   names=['chrom', 'start', 'end'])
        except FileNotFoundError:
            logging.warning(f"No real regions for {individual}, skipping.")
            continue

        real_df['length'] = real_df['end'] - real_df['start']
        # Extract per-chromosome real regions once
        per_chrom_real = {
            chrom: [(r.start, r.end) for r in IntervalTree(
                        Interval(row['start'], row['end']) 
                        for _, row in real_df[real_df['chrom']==chrom].iterrows()
                    )]
            for chrom in chromosomes
        }

        # Generate replicates
        for rep in range(1, num_replicates+1):
            bed_rows = []
            for chrom in chromosomes:
                chrom_regions = real_df[real_df['chrom'] == chrom]
                real_regions = [(row['start'], row['end']) for _, row in chrom_regions.iterrows()]
                # New interval tree for this replicate
                interval_tree = IntervalTree(Interval(s, e) for s, e in real_regions)
                for _, row in chrom_regions.iterrows():
                    region_length = row['length']
                    target_count = count_informative_snps(
                        snp_dict.get(chrom, []), row['start'], row['end']
                    )
                    try:
                        rs, re = generate_random_region(
                            chrom, region_length, target_count,
                            snp_dict.get(chrom, []), chrom_limits,
                            interval_tree, real_regions, snp_tol
                        )
                        bed_rows.append([chrom, rs, re])
                    except ValueError as e:
                        logging.warning(f"Rep {rep}: Unable to sample region for {chrom}:{row['start']}-{row['end']}: {e}")
                        continue

            # Write replicate BED
            out_dir = os.path.dirname(real_file)
            os.makedirs(out_dir, exist_ok=True)
            out_file = f"{individual}/{individual}.random.replicate{rep}_overlap.bed"
            pd.DataFrame(bed_rows, columns=['chrom', 'start', 'end']) \
                .to_csv(out_file, sep='\t', header=False, index=False)
            logging.info(f"Saved replicate {rep} for {individual} to {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate multiple sets of random regions with similar length and SNP count."
    )
    parser.add_argument("individuals_file", type=str,
                        help="File with individual names, one per line.")
    parser.add_argument('-n', '--num_replicates', type=int, default=100,
                        help="Number of random replicates to generate per individual.")
    args = parser.parse_args()
    main(args.individuals_file, args.num_replicates)

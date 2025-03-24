#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import logging
import csv
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from intervaltree import Interval, IntervalTree
import subprocess
import rpy2.robjects as robjects

# --------------------------------------------------------------------------------
# 1) Logging setup
# --------------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler()]
)

# --------------------------------------------------------------------------------
# 2) R code for density (example: exonic fraction)
# --------------------------------------------------------------------------------
r_script = r"""
library(GenomicRanges)
library(rtracklayer)

load_exons <- function(bed_path) {
    import(bed_path, format="bed")
}

calc_exonic_variant_fraction <- function(chrom, positions, bed_path) {
    exons <- load_exons(bed_path)
    exons_chrom <- exons[seqnames(exons) == chrom]

    if (length(positions) == 0) {
        return(0)
    }
    variants_gr <- GRanges(seqnames=chrom, IRanges(start=positions, end=positions))
    ov <- findOverlaps(variants_gr, exons_chrom)
    frac <- length(unique(queryHits(ov))) / length(positions)
    return(frac)
}
"""
robjects.r(r_script)
r_calc_exonic_fraction = robjects.globalenv['calc_exonic_variant_fraction']

def calculate_density_from_positions(chrom, positions, bed_path):
    """Wrapper to call R function to compute exonic fraction (or any density measure)."""
    if not positions:
        return 0.0
    positions_r = robjects.IntVector(positions)
    fraction_r = r_calc_exonic_fraction(chrom, positions_r, bed_path)
    return float(fraction_r[0])

# --------------------------------------------------------------------------------
# 3) bcftools extraction per chromosome (dictionary of VCF paths)
# --------------------------------------------------------------------------------
def get_variant_positions(vcf_dict, individual, chrom, start, end):
    """
    We have 22 separate VCF/BCF files, one per chromosome, stored in vcf_dict like:
       { "chr1": "/path/to/chr1.bcf", "chr2": "...", ..., "chr22": "..." }
    Use bcftools to extract variant positions in [start, end] for that chromosome+individual.
    """
    vcf_path = vcf_dict.get(chrom)
    if not vcf_path:
        logging.warning(f"No VCF path for {chrom}. Returning empty positions.")
        return []

    region_str = f"{chrom}:{start}-{end}"

    bcftools_view_cmd = [
        "bcftools", "view",
        "-s", individual,
        "-r", region_str,
        vcf_path,
        "-Ou"  # output uncompressed BCF
    ]
    bcftools_query_cmd = [
        "bcftools", "query",
        "-f", "%POS\\n"
    ]

    try:
        p1 = subprocess.Popen(bcftools_view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(bcftools_query_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        
        output, error = p2.communicate()
        if p2.returncode != 0:
            raise subprocess.CalledProcessError(p2.returncode, bcftools_view_cmd, output=error)
        
        lines = output.decode().strip().split('\n')
        positions = [int(x) for x in lines if x.strip()]
        return positions
    
    except FileNotFoundError:
        logging.error("bcftools not found on PATH.")
        return []
    except Exception as e:
        logging.error(f"bcftools error: {e}")
        return []

# --------------------------------------------------------------------------------
# 4) SNP counting (preloaded arrays)
# --------------------------------------------------------------------------------
def count_informative_snps_cached(snp_positions, start, end):
    if len(snp_positions) == 0:
        return 0
    left_idx = np.searchsorted(snp_positions, start, side='left')
    right_idx = np.searchsorted(snp_positions, end, side='right')
    return right_idx - left_idx

# --------------------------------------------------------------------------------
# 5) Random region generation
# --------------------------------------------------------------------------------
def generate_random_starts(range_length, min_start, max_start,
                           real_snp_count, snp_positions,
                           global_interval_tree,
                           snp_tol=0.05,
                           max_attempts_per_iteration=10000,
                           num_iterations=100):
    """
    Generate up to num_iterations random intervals that do not overlap any real region in global_interval_tree
    and have SNP count within ±5% of real_snp_count.
    """
    lower_bound = int(real_snp_count * (1 - snp_tol))
    upper_bound = int(real_snp_count * (1 + snp_tol))
    
    generated = []
    for _ in range(num_iterations):
        attempts = 0
        while attempts < max_attempts_per_iteration:
            rand_start = np.random.randint(min_start, max_start + 1)
            rand_end = rand_start + range_length

            # Check overlap with real intervals
            if global_interval_tree.overlap(rand_start, rand_end):
                attempts += 1
                continue

            # Check SNP count
            snp_count = count_informative_snps_cached(snp_positions, rand_start, rand_end)
            if lower_bound <= snp_count <= upper_bound:
                # Found valid region
                generated.append((rand_start, rand_end))
                break

            attempts += 1
        # If attempts == max_attempts_per_iteration => no region found => skip

    return generated

# --------------------------------------------------------------------------------
# 6) Worker function for a chunk of real regions
# --------------------------------------------------------------------------------
import logging

def process_region_chunk(chunk_rows, df_real_intervals, snp_dict, snp_tol, vcf_dict,
                         individual, exon_bed):
    """
    Process a subset of real regions: calculate densities for real and random regions.
    """
    process_id = os.getpid()
    logging.info(f"Process {process_id} started processing a chunk with {len(chunk_rows)} regions.")
    
    # Build the global IntervalTree containing all real regions
    global_tree = IntervalTree()
    for row in df_real_intervals.itertuples():
        global_tree.add(Interval(row.start, row.end))
    
    out_rows = []
    total_regions = len(chunk_rows)
    for idx, row in enumerate(chunk_rows.itertuples(), 1):
        chrom = row.chrom
        start = row.start
        end = row.end
        length = end - start
        real_region_str = f"{chrom}:{start}-{end}"
        
        # 1) Extract variant positions for the real region
        real_positions = get_variant_positions(vcf_dict, individual, chrom, start, end)
        real_density = calculate_density_from_positions(chrom, real_positions, exon_bed)
        
        # 2) Generate random regions with matched SNP counts
        snp_positions = snp_dict.get(chrom, np.array([]))
        real_snp_count = count_informative_snps_cached(snp_positions, start, end)
        min_start = row.first_pos
        max_start = row.last_pos - length
        
        random_intervals = generate_random_starts(
            range_length=length,
            min_start=min_start,
            max_start=max_start,
            real_snp_count=real_snp_count,
            snp_positions=snp_positions,
            global_interval_tree=global_tree,
            snp_tol=snp_tol,
            max_attempts_per_iteration=10000,
            num_iterations=100
        )
        
        # 3) Calculate density for each random region
        random_densities = []
        for (rstart, rend) in random_intervals:
            rand_positions = get_variant_positions(vcf_dict, individual, chrom, rstart, rend)
            rand_density = calculate_density_from_positions(chrom, rand_positions, exon_bed)
            random_densities.append(rand_density)
        
        # Fill up to 100 random densities if fewer were generated
        if len(random_densities) < 100:
            random_densities += [0.0] * (100 - len(random_densities))
        
        # 4) Compile the row for the CSV
        row_out = [real_region_str, real_density] + random_densities[:100]
        out_rows.append(row_out)
        
        # 5) Log progress every 10 regions
        if idx % 10 == 0 or idx == total_regions:
            logging.info(f"Process {process_id}: Processed {idx}/{total_regions} regions.")
    
    logging.info(f"Process {process_id} finished processing the chunk.")
    return out_rows

# --------------------------------------------------------------------------------
# 7) The parallel "process_individual" function
# --------------------------------------------------------------------------------
def process_individual_in_parallel(
    individual,
    vcf_dict,
    exon_bed,
    snp_dict,
    chrom_limits,
    snp_tolerance=0.05,
    n_cores=10
):
    """
    1) Load individual's .bed of real regions
    2) Attach first_pos, last_pos from chrom_limits
    3) Split into 'n_cores' chunks
    4) For each chunk, call 'process_region_chunk' in parallel
    5) Merge results in one CSV
    """
    import math
    
    input_bed = os.path.join(individual, f"{individual}.all.merged.tracts.bed")
    output_csv = os.path.join(individual, f"{individual}_variant_density.csv")

    try:
        df = pd.read_csv(input_bed, sep='\t', header=None, names=['chrom','start','end'])
    except Exception as e:
        logging.error(f"Error reading bed for {individual}: {e}")
        return

    df['original_length'] = df['end'] - df['start']
    df[['first_pos','last_pos']] = df['chrom'].map(
        lambda c: chrom_limits.get(c, [np.nan, np.nan])
    ).tolist()

    df = df.dropna(subset=['first_pos','last_pos'])

    n_rows = len(df)
    if n_rows == 0:
        logging.warning(f"No valid real regions found for {individual}.")
        return

    # We'll pass the entire DF as "df_real_intervals" so each chunk sees all intervals
    df_real_intervals = df.copy()

    # Divide the real regions into n_cores chunks
    chunk_size = int(np.ceil(n_rows / n_cores))
    chunks = []
    for i in range(0, n_rows, chunk_size):
        sub_df = df.iloc[i:i+chunk_size]
        chunks.append(sub_df)

    all_rows = []

    from concurrent.futures import ProcessPoolExecutor, as_completed
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = []
        for chunk_df in chunks:
            future = executor.submit(
                process_region_chunk,
                chunk_df,
                df_real_intervals,
                snp_dict,
                snp_tolerance,
                vcf_dict,
                individual,
                exon_bed
            )
            futures.append(future)

        for fut in as_completed(futures):
            try:
                result_rows = fut.result()
                all_rows.extend(result_rows)
            except Exception as e:
                logging.error(f"Chunk processing error: {e}")

    # Build header
    header = ["RealRegion", "RealDensity"] + [f"RandomDensity_{i+1}" for i in range(100)]

    # Save CSV
    os.makedirs(individual, exist_ok=True)
    try:
        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(all_rows)
        logging.info(f"Saved parallel results to {output_csv}")
    except Exception as e:
        logging.error(f"Error saving CSV: {e}")

# --------------------------------------------------------------------------------
# 8) Main
# --------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Parallel chunk-based script for real+random region densities, with multiple per-chrom VCFs.")
    parser.add_argument('individual', type=str, help="e.g. KSP092")
    parser.add_argument('--exon_bed', type=str, default="hg38.v47.exon.bed",
                        help="Path to BED file of exons or genes (for R overlap).")
    parser.add_argument('--chrom_limits', type=str, default="chrom_limits.txt",
                        help="File with chrom, first_pos, last_pos in tab-delimited format.")
    parser.add_argument('--snp_dir', type=str, default="SNP_files",
                        help="Directory containing files like chr1_informative_snps.txt, etc.")
    parser.add_argument('--cores', type=int, default=10, help="Number of parallel workers.")
    args = parser.parse_args()

    # Hard-coded 22 files, one per chromosome
    # Adapt paths as needed to your actual filenames.
    vcf_dict = {
        f"chr{i}": f"/lisc/scratch/admixlab/aigerim/ensembl-vep/sstar_san_{i}.vcf.gz" 
        for i in range(1, 23)
    }

    # Load chrom_limits
    try:
        df_limits = pd.read_csv(args.chrom_limits, sep='\t', header=None, names=['chrom','first_pos','last_pos'])
        chrom_limits = df_limits.set_index('chrom').T.to_dict('list')
    except Exception as e:
        logging.error(f"Failed to load chrom_limits: {args.chrom_limits} - {e}")
        return

    # Prepare a list of autosomes
    chromosomes = [f"chr{i}" for i in range(1, 23)]

    # Load informative SNP positions
    snp_dict = {}
    for c in chromosomes:
        cnum = c.replace('chr','')
        snp_file = os.path.join(args.snp_dir, f"chr{cnum}_informative_snps.txt")
        if not os.path.exists(snp_file):
            logging.warning(f"SNP file not found: {snp_file}, skipping {c}")
            snp_dict[c] = np.array([])
            continue
        try:
            arr = np.loadtxt(snp_file, usecols=1, dtype=int)
            snp_dict[c] = np.array(sorted(arr))
            logging.info(f"Loaded {len(arr)} SNP positions for {c}")
        except Exception as e:
            logging.error(f"Error loading {snp_file}: {e}")
            snp_dict[c] = np.array([])

    # Run the parallel chunk approach
    process_individual_in_parallel(
        individual=args.individual,
        vcf_dict=vcf_dict,
        exon_bed=args.exon_bed,
        snp_dict=snp_dict,
        chrom_limits=chrom_limits,
        snp_tolerance=0.05,
        n_cores=args.cores
    )

    logging.info(f"Done processing {args.individual}.")

if __name__ == "__main__":
    main()

import pandas as pd
import numpy as np
import subprocess
import multiprocessing
from functools import partial

def extract_sample_names(vcf_file):
    command = f"bcftools query -l {vcf_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    sample_names = result.stdout.split()
    return sample_names

def count_diff_with_groups(bed_file, vcf_file, group1_ids, group2_ids, individual):
    # Load BED file with regions
    regions = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
    
    # Initialize an empty DataFrame for aggregated SNP data
    aggregated_data = []

    # Process each region specified in the BED file
    for index, row in regions.iterrows():
        region = f"{row['chrom']}:{row['start']}-{row['end']}"
        command = f"bcftools view -m2 -M2 -v snps -r {region} {vcf_file} | bcftools query -f '%CHROM %POS %REF %ALT [ %GT ]\\n'"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.stdout:
            testsnps = result.stdout.splitlines()
            data = [line.split() for line in testsnps]
            aggregated_data.extend(data)

    if not aggregated_data:
        return "No SNP data found in the specified regions."

    # Extract sample names from the VCF file
    sample_names = extract_sample_names(vcf_file)

    # Convert aggregated SNP data into a DataFrame
    aut_oneid = pd.DataFrame(aggregated_data)

    # Set DataFrame column indices for genotypes, using sample names from the VCF file
    aut_oneid.columns = ['CHROM', 'POS', 'REF', 'ALT'] + sample_names

    # Print DataFrame shape before processing
    print(f"Shape before processing: {aut_oneid.shape}")

    # Convert genotype entries, treating "./." as a special case
    def convert_genotypes(gt):
        if gt == "0/0":
            return 0
        elif gt in ["0/1", "1/0"]:
            return np.random.choice([0, 1])  # Randomly assign 0 or 1 for heterozygotes
        elif gt == "1/1":
            return 1
        elif gt == "./.":
            return np.nan  # Treat missing data as -1
        return np.nan  # Handle unexpected genotypes

    for col in aut_oneid.columns[4:]:  # Skip chromosomal position data
        aut_oneid[col] = aut_oneid[col].apply(convert_genotypes)
    
    # Print DataFrame shape after converting genotypes
    print(f"Shape after converting genotypes: {aut_oneid.shape}")

    # Remove rows with any missing data (./.)
    aut_oneid.dropna(inplace=True)

    # Print DataFrame shape after removing "./." data
    print(f"Shape after removing missing data: {aut_oneid.shape}")

    # Extract the relevant genotype data using sample names
    group1 = aut_oneid[group1_ids]
    group2 = aut_oneid[group2_ids]

    # Function to calculate differences between an individual and all individuals in a group
    def calculate_individual_differences(individual, group1, group2):
        # Check if the individual exists in group1
        if individual not in group1.columns:
            raise ValueError(f"Individual {individual} not found in group1")

        differences = {}

        # Calculate differences with individuals in group1
        for col in group1.columns:
            if col != individual:
                diff_sum = (group1[individual] != group1[col]).sum()
                differences[col] = diff_sum

        # Calculate differences with individuals in group2
        for col in group2.columns:
            diff_sum = (group1[individual] != group2[col]).sum()
            differences[col] = diff_sum

        return differences

    # Calculate total pairwise differences for the given individual
    result = calculate_individual_differences(individual, group1, group2)

    return result

def process_individual(individual, group1_ids, group2_ids):
    bed_file = f"pairwisediff/{individual}_random_overlap.bed"
    total_diffs = {ind: 0 for ind in group1_ids + group2_ids if ind != individual}

    for chrom in range(1, 23):
        vcf_file = f"25KS.48RHG.74comp.HCBP.{chrom}.recalSNP99.9.recalINDEL99.0.vcf.gz"
        result = count_diff_with_groups(bed_file, vcf_file, group1_ids, group2_ids, individual)
        if isinstance(result, dict):
            for ind, diff_sum in result.items():
                total_diffs[ind] += diff_sum
    return individual, total_diffs

if __name__ == "__main__":
    output_directory = "pairwisediff"
    group1_ids = ["KSP062", "KSP063", "KSP065", "KSP067", "KSP069", "KSP092", "KSP096", "KSP103", "KSP105", "KSP106", "KSP111", "KSP116", "KSP124", "KSP134", "KSP137", "KSP139", "KSP140", "KSP146", "KSP150", "KSP152", "KSP154", "KSP155", "KSP224", "KSP225", "KSP228"] # target
    group2_ids = ["HG02568", "HG02922", "HG03052", "HGDP_HGDP00927", "HGDP_HGDP01284", "SGDP_LP6005441-DNA_E07", "SGDP_LP6005441-DNA_F07", "SGDP_LP6005442-DNA_A02", "SGDP_LP6005442-DNA_A10", "SGDP_LP6005442-DNA_B02", "SGDP_LP6005442-DNA_B10", "SGDP_LP6005442-DNA_G10", "SGDP_LP6005442-DNA_G11", "SGDP_LP6005442-DNA_H10", "SGDP_SS6004475", "SGDP_SS6004470", "HGDP_DNK02", "NA19017", "PNP010", "PNP011", "PNP012", "PNP013", "PNP014", "PNP030", "PNP031"]

    individuals = ["HGDP_HGDP01029", "KSP062", "KSP063", "KSP065", "KSP067", "KSP069", "KSP092", "KSP096", "KSP103", "KSP105", "KSP106", "KSP111", "KSP116", "KSP124", "KSP134", "KSP137", "KSP139", "KSP140", "KSP146", "KSP150", "KSP152", "KSP154", "KSP155", "KSP224", "KSP225", "KSP228"]

    final_results_df = pd.DataFrame()

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(partial(process_individual, group1_ids=group1_ids, group2_ids=group2_ids), individuals)

    for individual, total_diffs in results:
        total_diff_df = pd.DataFrame.from_dict(total_diffs, orient='index', columns=[individual])
        if final_results_df.empty:
            final_results_df = total_diff_df
        else:
            final_results_df = final_results_df.join(total_diff_df, how='outer')

    final_diff_csv_path = f"{output_directory}/random_diff_sums_target_individuals.csv"
    final_results_df.to_csv(final_diff_csv_path, index_label='Individual')
    print(f"Total sum of differences for all individuals saved to: {final_diff_csv_path}")

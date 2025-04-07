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

def count_diff_with_group2(bed_file, vcf_file, group2_ids, individual2, max_rows=287605):
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

    # Convert genotype entries, treating "./." as NaN
    def convert_genotypes(gt):
        if gt == "0/0":
            return 0
        elif gt in ["0/1", "1/0"]:
            return np.random.choice([0, 1])  # Randomly assign 0 or 1 for heterozygotes
        elif gt == "1/1":
            return 1
        elif gt == "./.":
            return np.nan  # Treat missing data as NaN
        return np.nan  # Handle unexpected genotypes

    for col in aut_oneid.columns[4:]:  # Skip chromosomal position data
        aut_oneid[col] = aut_oneid[col].apply(convert_genotypes)

    # Print DataFrame shape after converting genotypes
    print(f"Shape after converting genotypes: {aut_oneid.shape}")

    # Remove rows with any missing data (./. handled as NaN)
    aut_oneid.dropna(inplace=True)

    # Print DataFrame shape after removing rows with NaN values
    print(f"Shape after removing rows with NaN: {aut_oneid.shape}")

    # Restrict DataFrame to max_rows
    if len(aut_oneid) > max_rows:
        aut_oneid = aut_oneid.head(max_rows)

    # Print DataFrame shape after restricting rows
    print(f"Shape after restricting to {max_rows} rows: {aut_oneid.shape}")

    # Extract the relevant genotype data using sample names
    group2 = aut_oneid[group2_ids]

    # Function to calculate differences between an individual and all individuals in group2
    def calculate_individual_differences(individual2, group2):
        # Check if the individual exists in group2
        if individual2 not in group2.columns:
            raise ValueError(f"Individual {individual2} not found in group2")

        differences = {}

        # Calculate differences with individuals in group2
        for col in group2.columns:
            if col != individual2:
                diff_sum = (group2[individual2] != group2[col]).sum()
                differences[col] = diff_sum

        return differences

    # Calculate total pairwise differences for the given individual
    result = calculate_individual_differences(individual2, group2)

    return result

def process_individual(individual2, bed_files, group2_ids, max_rows=287605):
    # Get the BED file for the current individual
    bed_file = bed_files.get(individual2)
    if not bed_file:
        print(f"No BED file found for {individual2}, skipping.")
        return individual2, {}

    # Initialize totals for all chromosomes for the current individual
    total_diffs = {ind: 0 for ind in group2_ids if ind != individual2}

    # Loop through chromosome numbers 1 to 22
    for chrom in range(1, 23):
        vcf_file = f"25KS.48RHG.74comp.HCBP.{chrom}.recalSNP99.9.recalINDEL99.0.vcf.gz"
        result = count_diff_with_group2(bed_file, vcf_file, group2_ids, individual2, max_rows)
        if isinstance(result, str):  # Check if result is an error message
            print(f"Chromosome {chrom} result for {individual2}: {result}")
            continue
        print(f"Chromosome {chrom} result for {individual2}: {result}")
        for ind, diff_sum in result.items():
            total_diffs[ind] += diff_sum
            print(f"Total differences for {individual2} with {ind}: {total_diffs[ind]}")

    return individual2, total_diffs

if __name__ == "__main__":
    output_directory = "pairwisediff"
    bed_files = {
        "HG02922": "pairwisediff/KSP062_random_overlap.bed",
        "HG03052": "pairwisediff/KSP063_random_overlap.bed",
        "HGDP_HGDP00927": "pairwisediff/KSP065_random_overlap.bed",
        "HGDP_HGDP01284": "pairwisediff/KSP067_random_overlap.bed",
        "SGDP_LP6005441-DNA_E07": "pairwisediff/KSP069_random_overlap.bed",
        "SGDP_LP6005441-DNA_F07": "pairwisediff/KSP092_random_overlap.bed",
        "SGDP_LP6005442-DNA_A02": "pairwisediff/KSP096_random_overlap.bed",
        "SGDP_LP6005442-DNA_A10": "pairwisediff/KSP103_random_overlap.bed",
        "SGDP_LP6005442-DNA_B02": "pairwisediff/KSP105_random_overlap.bed",
        "SGDP_LP6005442-DNA_B10": "pairwisediff/KSP106_random_overlap.bed",
        "SGDP_LP6005442-DNA_G10": "pairwisediff/KSP111_random_overlap.bed",
        "SGDP_LP6005442-DNA_G11": "pairwisediff/KSP116_random_overlap.bed",
        "SGDP_LP6005442-DNA_H10": "pairwisediff/KSP124_random_overlap.bed",
        "SGDP_SS6004475": "pairwisediff/KSP134_random_overlap.bed",
        "SGDP_SS6004470": "pairwisediff/KSP137_random_overlap.bed",
        "HGDP_DNK02": "pairwisediff/KSP139_random_overlap.bed",
        "NA19017": "pairwisediff/KSP140_random_overlap.bed",
        "PNP010": "pairwisediff/KSP146_random_overlap.bed",
        "PNP011": "pairwisediff/KSP150_random_overlap.bed",
        "PNP012": "pairwisediff/KSP152_random_overlap.bed",
        "PNP013": "pairwisediff/KSP154_random_overlap.bed",
        "PNP014": "pairwisediff/KSP155_random_overlap.bed",
        "PNP030": "pairwisediff/KSP224_random_overlap.bed",
        "PNP031": "pairwisediff/KSP225_random_overlap.bed",
        "PNP032": "pairwisediff/KSP228_random_overlap.bed"
    }
    group2_ids = ["HG02568", "HG02922", "HG03052", "HGDP_HGDP00927", "HGDP_HGDP01284", "SGDP_LP6005441-DNA_E07", "SGDP_LP6005441-DNA_F07", "SGDP_LP6005442-DNA_A02", "SGDP_LP6005442-DNA_A10", "SGDP_LP6005442-DNA_B02", "SGDP_LP6005442-DNA_B10", "SGDP_LP6005442-DNA_G10", "SGDP_LP6005442-DNA_G11", "SGDP_LP6005442-DNA_H10", "SGDP_SS6004475", "SGDP_SS6004470", "HGDP_DNK02", "NA19017", "PNP010", "PNP011", "PNP012", "PNP013", "PNP014", "PNP030", "PNP031", "PNP032"]

   # List of individuals to process
    individuals = ["HG02568", "HG02922", "HG03052", "HGDP_HGDP00927", "HGDP_HGDP01284", "SGDP_LP6005441-DNA_E07", "SGDP_LP6005441-DNA_F07", "SGDP_LP6005442-DNA_A02", "SGDP_LP6005442-DNA_A10", "SGDP_LP6005442-DNA_B02", "SGDP_LP6005442-DNA_B10", "SGDP_LP6005442-DNA_G10", "SGDP_LP6005442-DNA_G11", "SGDP_LP6005442-DNA_H10", "SGDP_SS6004475", "SGDP_SS6004470", "HGDP_DNK02", "NA19017", "PNP010", "PNP011", "PNP012", "PNP013", "PNP014", "PNP030", "PNP031", "PNP032"]

    final_results_df = pd.DataFrame()

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(partial(process_individual, bed_files=bed_files, group2_ids=group2_ids), individuals)

    for individual2, total_diffs in results:
        if not total_diffs:
            continue
        total_diff_df = pd.DataFrame.from_dict(total_diffs, orient='index', columns=[individual2])
        if final_results_df.empty:
            final_results_df = total_diff_df
        else:
            final_results_df = final_results_df.join(total_diff_df, how='outer')

    final_diff_csv_path = f"{output_directory}/min_random_diff_sums_all_ref_individuals.csv"
    final_results_df.to_csv(final_diff_csv_path, index_label='Individual')
    print(f"Total sum of differences for all individuals saved to: {final_diff_csv_path}")

    # Print the final results
    print("Final total sum of differences for all individuals:")
    print(final_results_df)


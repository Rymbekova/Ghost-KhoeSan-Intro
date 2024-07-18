import numpy as np
import pandas as pd
import pyranges as pr

def read_genomic_data(filepath):
    """Read genomic data from a file."""
    df = pd.read_csv(filepath, sep="\t", header=None, names=['Chromosome', 'Start', 'End'])
    print("Data read from file:", filepath)
    #print(df.head())
    return pr.PyRanges(df)

def read_chromosome_lengths(filepath):
    """Read chromosome lengths from a file."""
    chrom_lengths = pd.read_csv(filepath, sep="\t", header=None, names=['Chromosome', 'Length'])
    chrom_dict = chrom_lengths.set_index('Chromosome')['Length'].to_dict()
    print("Chromosome lengths read:")
    #print(chrom_dict)
    return chrom_dict

def generate_random_genomic_data_by_chromosome(chrom_lengths, original_data):
    """Generate random genomic regions based on original data's chromosome distribution."""
    data_frames = []
    for chrom, group in original_data.df.groupby('Chromosome'):
        length = chrom_lengths.get(chrom, 0)
        count = len(group)
        lengths = group['End'] - group['Start']
        starts = np.random.randint(0, length - lengths.max(), size=count)
        ends = starts + lengths
        df = pd.DataFrame({'Chromosome': [chrom]*count, 'Start': starts, 'End': ends})
        data_frames.append(df)
        #print(f"Generated random data for {chrom}:")
        #print(df.head())

    combined_data = pd.concat(data_frames)
    return pr.PyRanges(combined_data)

def calculate_overlaps(data1, data2):
    """Calculate and return the overlap as a percentage of the base pairs for both data sets."""
    overlaps = data1.intersect(data2)
    total_overlap_length = overlaps.lengths().sum()
    total_length_data1 = data1.lengths().sum()
    total_length_data2 = data2.lengths().sum()
    overlap_percent_data1 = (total_overlap_length / total_length_data1) * 100 if total_length_data1 else 0
    overlap_percent_data2 = (total_overlap_length / total_length_data2) * 100 if total_length_data2 else 0
    print(f"Overlap percent for this iteration (sstar): {overlap_percent_data1}%")
    print(f"Overlap percent for this iteration (skov): {overlap_percent_data2}%")
    return overlap_percent_data1, overlap_percent_data2

# Main execution block
if __name__ == "__main__":
    filepath1 = 'KSP228/KSP228.all.merged.tracts.bed'
    filepath2 = '../Introgression-detection/KSP228/KSP228.all.merged.90.bed'
    chrom_lengths_path = 'chrom_sizes.txt'
    iterations = 100
    results_data1 = []
    results_data2 = []

    print("Starting process...")
    data1 = read_genomic_data(filepath1)
    data2 = read_genomic_data(filepath2)
    chrom_lengths = read_chromosome_lengths(chrom_lengths_path)

    for i in range(iterations):
        print(f"--- Iteration {i+1} ---")
        random_data1 = generate_random_genomic_data_by_chromosome(chrom_lengths, data1)
        random_data2 = generate_random_genomic_data_by_chromosome(chrom_lengths, data2)
        
        overlap_percent_data1, overlap_percent_data2 = calculate_overlaps(random_data1, random_data2)
        results_data1.append(overlap_percent_data1)
        results_data2.append(overlap_percent_data2)

    average_overlap_data1 = np.mean(results_data1)
    average_overlap_data2 = np.mean(results_data2)
    print("Average Overlap Across 100 Iterations (sstar):", average_overlap_data1)
    print("Average Overlap Across 100 Iterations (skov):", average_overlap_data2)

import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree

# Load data without headers
ranges_df = pd.read_csv('HGDP_HGDP01029/HGDP_HGDP01029_overlap.bed', sep='\t', header=None, names=['chrom', 'start', 'end'])
limits_df = pd.read_csv('chrom_limits.txt', sep='\t', header=None, names=['chrom', 'first_pos', 'last_pos'])

# Calculate length of original ranges
ranges_df['original_length'] = ranges_df['end'] - ranges_df['start']

# Merge limits into ranges to ensure matching chromosomes are together
ranges_df = ranges_df.merge(limits_df, on='chrom', how='left')

# Check for missing limits
missing_limits = ranges_df[ranges_df['first_pos'].isnull()]
if not missing_limits.empty:
    print("Warning: Some ranges have no corresponding chromosomal limits and will be excluded:")
    print(missing_limits)
    ranges_df = ranges_df.dropna(subset=['first_pos'])

# Create a dictionary of interval trees for each chromosome
chrom_interval_trees = {chrom: IntervalTree() for chrom in ranges_df['chrom'].unique()}

# Populate the interval trees with existing intervals
for chrom, start, end in zip(ranges_df['chrom'], ranges_df['start'], ranges_df['end']):
    chrom_interval_trees[chrom].add(Interval(start, end))

# Function to generate a random start position within chromosomal limits and avoiding overlaps
def generate_random_start(row, chrom_interval_trees):
    range_length = row['end'] - row['start']
    min_start = row['first_pos']
    max_start = row['last_pos'] - range_length
    max_attempts = 1000
    attempts = 0
    chrom = row['chrom']
    
    while attempts < max_attempts:
        random_start = np.random.randint(min_start, max_start + 1)
        random_end = random_start + range_length
        # Check for overlap with any existing intervals in the same chromosome using the interval tree
        if not chrom_interval_trees[chrom].overlap(random_start, random_end):
            # Add new interval to the tree to ensure no future overlaps
            chrom_interval_trees[chrom].add(Interval(random_start, random_end))
            return random_start
        attempts += 1
    
    # Return NaN if no valid start found after max_attempts
    return np.nan

# Apply function and calculate new end, while checking for overlaps
ranges_df['new_start'] = ranges_df.apply(generate_random_start, axis=1, chrom_interval_trees=chrom_interval_trees)
ranges_df = ranges_df.dropna(subset=['new_start'])
ranges_df['new_end'] = ranges_df['new_start'] + (ranges_df['end'] - ranges_df['start'])

# Create output DataFrame with new ranges
output_df = ranges_df[['chrom', 'new_start', 'new_end']]

# Sort the output DataFrame by chromosome and start position
output_df = output_df.sort_values(by=['chrom', 'new_start'])

# Save to file
output_df.to_csv('HGDP_HGDP01029/HGDP_HGDP01029_random_overlap.bed', sep='\t', index=False, header=False)

# Verification steps
# Load the generated file
output_df = pd.read_csv('HGDP_HGDP01029/HGDP_HGDP01029_random_overlap.bed', sep='\t', header=None, names=['chrom', 'start', 'end'])

# Calculate the length of each new range
output_df['new_length'] = output_df['end'] - output_df['start']

# Calculate the total length of all original ranges
total_original_length = ranges_df['original_length'].sum()

# Calculate the total length of all new ranges
total_new_length = output_df['new_length'].sum()

# Print the shape of the output DataFrame
print("\nShape of the output DataFrame:", output_df.shape)

# Print the head of the output DataFrame with the length of each new range and the total lengths
print("Head of the output DataFrame with lengths of each new range:")
print(output_df.head())

print("\nTotal length of all original ranges:", total_original_length)
print("Total length of all new ranges:", total_new_length)

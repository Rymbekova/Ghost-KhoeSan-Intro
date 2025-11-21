def calculate_total_length(file_path):
    """Calculate the total length of genomic fragments in a BED file."""
    total_length = 0
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            start = int(parts[1])
            end = int(parts[2])
            fragment_length = end - start
            total_length += fragment_length
    return total_length

# List of base paths or patterns to generate file names for each kind
base_paths = [
    'KSP228/KSP228.{}.ref50.tgt1.merged.tracts.bed',
    '../Introgression-detection/KSP228/KSP228.{}.ref50.tgt1.recode.diploid.90.bed'
]

# Dictionary to keep total lengths for two kinds of files: sstar and skov
total_lengths = {1: 0, 2: 0}
for i in range(1, 23):  # Assuming chromosome numbers or identifiers from 1 to 22
    for j, base_path in enumerate(base_paths, start=1):
        file_path = base_path.format(i)
        total_length = calculate_total_length(file_path)
        total_lengths[j] += total_length
        file_type = 'sstar' if j == 1 else 'skov'
        #print(f"Total length of fragments in {file_type} file {file_path}:", total_length)

# Output the total lengths summed across all files for each kind
print("Total length of all sstar files:", total_lengths[1])
print("Total length of all skov files:", total_lengths[2])

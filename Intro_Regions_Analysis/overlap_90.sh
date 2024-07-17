#!/bin/bash

cd ../Introgression-detection/KSP062/

# Loop over each .90.bed file
for f in *.diploid.95.bed; do
    # Extract the chromosome number (e.g., 2)
    chromosome=$(echo "$f" | awk -F'.' '{print $2}')

    # Construct the name of the corresponding merged tracts file
    merged_file="../../San21/KSP062/KSP062.$chromosome.ref50.tgt1.merged.tracts.bed"

    # Calculate the overlap using bedtools
    overlap=$(bedtools intersect -a "$f" -b "$merged_file" -wo | awk '{sum += $NF} END {print sum}')

    # Echo the chromosome number and the calculated overlap
    echo "Chromosome: $chromosome, Overlap: $overlap"
done


# Initialize a variable to keep the total overlap sum
total_overlap=0

# Loop over each .90.bed file
for f in *.diploid.95.bed; do
    # Extract the chromosome number (e.g., 2 from HGDP_HGDP01029.2...)
    chromosome=$(echo "$f" | awk -F'.' '{print $2}')

    # Construct the name of the corresponding merged tracts file
    merged_file="../../San21/KSP062/KSP062.$chromosome.ref50.tgt1.merged.tracts.bed"

    # Calculate the overlap using bedtools
    overlap=$(bedtools intersect -a "$f" -b "$merged_file" -wo | awk '{sum += $NF} END {print sum}')

    # Add the overlap for this pair to the total overlap
    total_overlap=$(($total_overlap + overlap))
done

# Print the total overlap across all chromosomes
echo "KSP062 overlap all chroms at 90%: $total_overlap"

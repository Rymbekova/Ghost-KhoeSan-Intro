#!/bin/bash

cd ../Introgression-detection/HGDP_HGDP01029/

output_dir="../../San21/HGDP_HGDP01029"

# Loop over each .90.bed file
for f in *.diploid.90.bed; do
    # Extract the chromosome number
    chromosome=$(echo "$f" | awk -F'.' '{print $2}')

    # Construct the name of the corresponding merged tracts file
    merged_file="../../San21/HGDP_HGDP01029/HGDP_HGDP01029.$chromosome.ref50.tgt1.merged.tracts.bed"

    # Use bedtools to find the overlaps and write them to a file
    output_file="${output_dir}/chromosome_${chromosome}_overlap.bed"
    bedtools intersect -a "$f" -b "$merged_file" -wo | \
    awk 'BEGIN {OFS="\t"} {start = $2 > $5 ? $2 : $5; overlap_size = $NF; print $1, start, start + overlap_size}' > "$output_file"

    # Echo the chromosome number and the file where overlaps were written
    echo "Chromosome: $chromosome, Overlaps written to: $output_file"
done

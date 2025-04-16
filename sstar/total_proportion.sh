#!/bin/bash

# Check if the input file exists
input_file="all.San.ind"

module load bedtools

if [ ! -f "$input_file" ]; then
    echo "Input file $input_file not found."
    exit 1
fi

# Loop through sample names from the input file
while IFS= read -r sample_name; do
    # Initialize variables to store merged intervals
    score_merged_file="${sample_name}_score_merged.bed"
    tracts_merged_file="${sample_name}_tracts_merged.bed"
    
    # Merge intervals from the score and tracts BED files for all chromosomes
    for chr in {1..22}; do
	awk '$5 != "NA"' "$sample_name/$sample_name.$chr.ref50.tgt1.score.results" > "$sample_name/$sample_name.$chr.ref50.tgt1.filtered.score.results"
        bedtools merge -i "$sample_name/$sample_name.$chr.ref50.tgt1.filtered.score.results" >> "$score_merged_file"
        bedtools merge -i "$sample_name/$sample_name.$chr.ref50.tgt1.tracts.bed" >> "$tracts_merged_file"
    done

    # Calculate the length of the merged score and tracts intervals
    score_length=$(awk '{sum += ($3 - $2)} END {print sum}' "$score_merged_file")
    echo "$score_length"
    tracts_length=$(awk '{sum += ($3 - $2)} END {print sum}' "$tracts_merged_file")
    echo "$tracts_length"

    # Calculate the proportion with higher precision
    if [ "$score_length" -ne 0 ]; then
        proportion=$(bc <<< "scale=4; ($tracts_length / $score_length) * 100")
    else
        proportion="N/A"  # Avoid division by zero
    fi

    # Print the result with higher precision
    echo "Proportion for $sample_name: $proportion"

    # Optionally, you can remove temporary merged BED files if needed
    rm "$score_merged_file"
    rm "$tracts_merged_file"
    rm "$sample_name/$sample_name.$chr.ref50.tgt1.filtered.score.results"

done < "$input_file"


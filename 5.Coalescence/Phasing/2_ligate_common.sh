#!/bin/bash
#SBATCH --job-name=ligate
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=3:00:00
#SBATCH --output=./pipeline_ligate.out
#SBATCH --error=./pipeline_ligate.err

module load glimpse/2.0.0
module load htslib
module load libdeflate/1.24

ligate_exe="../ligate_static"

input_dir="chunks_common"  
output_dir="ligate_common"
thread=2

mkdir -p "$output_dir"

# Loop through the chromosomes
for chr in {1..22}; do
    echo "Processing chr${chr}..."
    
    # Generate the ordered list of chunks for the current chromosome
    chunk_list="$output_dir/list_chr${chr}.txt"
    find "$input_dir" -name "Archaic_African_chr${chr}.chunk_*.shapeit5_common.bcf" | sort -V > "$chunk_list"
    
    # Define output file for the current chromosome
    output_file="$output_dir/Archaic_African_chr${chr}.ligated.bcf"
    
    # Run ligate command
    "$ligate_exe" --input "$chunk_list" \
                  --output "$output_file" \
                  --thread $thread
done

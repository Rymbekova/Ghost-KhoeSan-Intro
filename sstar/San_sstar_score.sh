#!/bin/bash
#SBATCH --job-name=sstar_score
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=36:00:00

module load conda
conda activate sstar-analysis

# Check if the input file exists
input_file="all.San.ind"

if [ ! -f "$input_file" ]; then
    echo "Input file $input_file not found."
    exit 1
fi

# Read sample names from the input file
while IFS= read -r sample_name; do
    # Navigate to the folder
    cd "$sample_name"
    
    # Run the sstar score command
    sstar score --vcf "${sample_name}.ref50.tgt1.recode.vcf" --ref ref.ind --tgt tgt.ind --output "${sample_name}.ref50.tgt1.score.results"
    
    echo "sstar score command for $sample_name executed successfully."
    
    # Return to the parent directory
    cd ..
done < "$input_file"

echo "sstar score commands executed successfully for all sample name folders."


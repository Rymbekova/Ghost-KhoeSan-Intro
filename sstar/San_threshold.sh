#!/bin/bash
#SBATCH --job-name=San_threshold
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=4:00:00

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

    for chr in {1..22}; do 
        sstar threshold --score "$sample_name/$sample_name.$chr.ref50.tgt1.score.results" --sim-data /lisc/scratch/admixlab/aigerim/sstar-analysis-main/results/inference/sstar/AfricanNullSan/quantiles/quantile.summary.txt  --quantile 0.99 --output "$sample_name/$sample_name.$chr.ref50.tgt1.threshold.results"
        echo "sstar threshold command for $sample_name ($chr) executed successfully."
        sstar tract --threshold "$sample_name/$sample_name.$chr.ref50.tgt1.threshold.results" --output-prefix "$sample_name/$sample_name.$chr.ref50.tgt1.tracts"
        echo "sstar tract command for $sample_name ($chr) executed successfully."
    done
done < "$input_file"

echo "sstar commands executed successfully for all sample names."

#!/bin/bash
#SBATCH --job-name=ref_random
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=./pipeline_ref_random.out
#SBATCH --error=./pipeline_ref_random.err

module load conda
conda activate sstar-analysis
module load bcftools

python random_ref_sum.py


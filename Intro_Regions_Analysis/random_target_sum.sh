#!/bin/bash
#SBATCH --job-name=target_random
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=./pipeline_target_random.out
#SBATCH --error=./pipeline_target_random.err

module load conda
conda activate sstar-analysis
module load bcftools

python random_target_sum.py


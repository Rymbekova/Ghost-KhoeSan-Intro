#!/bin/bash
#SBATCH --job-name=target
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=./pipeline_target.out
#SBATCH --error=./pipeline_target.err

module load conda
conda activate sstar-analysis
module load bcftools

python target_sum.py


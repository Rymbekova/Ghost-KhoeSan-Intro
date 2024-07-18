#!/bin/bash
#SBATCH --job-name=ref
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=./pipeline_ref.out
#SBATCH --error=./pipeline_ref.err

module load conda
conda activate sstar-analysis
module load bcftools

python ref_sum.py


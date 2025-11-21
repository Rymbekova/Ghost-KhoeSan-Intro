#!/bin/bash
#SBATCH --job-name=quantile
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=24:00:00

module load conda
conda activate sstar-analysis

sstar quantile --model config/simulation/models/null_african.yaml --ms-dir ext/msdir --N0 1000 --nsamp 102 --nreps 20000 --ref-index 5 --ref-size 100 --tgt-index 7 --tgt-size 2 --mut-rate 1.29e-8 --rec-rate 1e-8 --seq-len 40000 --snp-num-range 25 705 5 --output-dir results/inference/sstar/AfricanNullSan/quantiles --thread 16


#!/bin/bash
#SBATCH --job-name=concat_rare
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=4:00:00
#SBATCH --array=2-22
#SBATCH --output=logs/concat_rare_%A_%a.out
#SBATCH --error=logs/concat_rare_%A_%a.err

module load bcftools

CHR=${SLURM_ARRAY_TASK_ID}
threads=4

echo "[chr${CHR}] Starting concatenation..."

ODIR="final_concat"
mkdir -p "$ODIR"

find chunks_rare/ -name "Archaic_African_chr${CHR}.chunk_*.shapeit5_rare.bcf" | sort -V > "${ODIR}/concat_list_chr${CHR}.txt"

OUT="${ODIR}/Archaic_African_chr${CHR}.full.shapeit5_rare.bcf"

bcftools concat -n -f "${ODIR}/concat_list_chr${CHR}.txt" -o "$OUT" -O b
bcftools index "$OUT" --threads "$threads"

rm "${ODIR}/concat_list_chr${CHR}.txt"

echo "[chr${CHR}] Done."

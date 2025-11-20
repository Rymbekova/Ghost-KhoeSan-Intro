#!/bin/bash
#SBATCH --job-name=phase_chr_batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --array=1-22
#SBATCH --output=logs/phase_chr_batch_%A_%a.out
#SBATCH --error=logs/phase_chr_batch_%A_%a.err

module load glimpse/2.0.0
module load htslib/1.20

phase_common_exe="/lisc/scratch/admixlab/aigerim/shapeit5/phase_common/bin/phase_common"
threads=4
maf=0.001
PREFIX="Archaic_African"

CHR=${SLURM_ARRAY_TASK_ID}
CHUNKS="chunks/phase_common/chunks_chr${CHR}.txt"
BCF="bcf/Archaic_African_chr${CHR}.bcf"
MAP="/lisc/scratch/admixlab/aigerim/glimpse/maps/genetic_maps.b38/chr${CHR}.b38.with_chr.gmap.gz"

while read -r LINE; do
    CHUNK_NBR=$(echo "$LINE" | awk '{ print $1 }')
    REGION=$(echo "$LINE" | awk '{ print $3 }')
    OUTPUT_FILE="chunks_common/${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf"

    "$phase_common_exe" \
        --input "$BCF" \
        --map "$MAP" \
        --output "$OUTPUT_FILE" \
        --thread "$threads" \
        --filter-maf "$maf" \
        --region "$REGION"

done < "$CHUNKS"

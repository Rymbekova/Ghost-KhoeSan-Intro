#!/bin/bash
#SBATCH --job-name=phase_rare_chr_batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=3:00:00
#SBATCH --array=7-22
#SBATCH --output=logs/phase_rare_chr_batch_%A_%a.out
#SBATCH --error=logs/phase_rare_chr_batch_%A_%a.err

module load glimpse/2.0.0
module load htslib/1.20

phase_rare_exe="../phase_rare_static"
threads=4

CHR=${SLURM_ARRAY_TASK_ID}
CHUNKS="../African/chunks/phase_rare/chunks_chr${CHR}.txt"
BCF="Archaic_African.chr${CHR}.bcf"
MAP="/lisc/scratch/admixlab/aigerim/glimpse/maps/genetic_maps.b38/chr${CHR}.b38.with_chr.gmap.gz"
SCAF="ligate_common/Archaic_African_chr${CHR}.ligated.bcf"
ODIR="chunks_rare"
mkdir -p "$ODIR"

while read -r LINE; do
    CHUNK_NBR=$(echo "$LINE" | awk '{ print $1 }')
    SCAFFOLD_REG=$(echo "$LINE" | awk '{ print $3 }')
    INPUT_REG=$(echo "$LINE" | awk '{ print $4 }')
    OUTPUT_FILE="${ODIR}/Archaic_African_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.bcf"

    "$phase_rare_exe" \
        --input "$BCF" \
        --map "$MAP" \
        --output "$OUTPUT_FILE" \
        --thread "$threads" \
        --scaffold "$SCAF" \
        --scaffold-region "$SCAFFOLD_REG" \
        --input-region "$INPUT_REG"

done < "$CHUNKS"

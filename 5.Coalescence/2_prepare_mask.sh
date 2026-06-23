#!/bin/bash
#SBATCH --job-name=prepare_relate_inputs
#SBATCH --mem=300G
#SBATCH --time=3:00:00
#SBATCH --array=1-22
#SBATCH --output=prep_chr%a.log

PREP_SCRIPT="/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh"

chr=${SLURM_ARRAY_TASK_ID}

# Input haps & sample 
haps_file="haps/Archaic_African_chr${chr}.haps"
sample_file="sample/Archaic_African_${chr}.sample"

# Mask & ancestor files
mask_file="/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/African/Relate_input_files/GRCh38/20160622_genome_mask_GRCh38/PilotMask/20160622.chr${chr}.pilot_mask.fasta.gz"
ancestor_file="/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/African/Relate_input_files/GRCh38/human_ancestor_GRCh38/homo_sapiens_ancestor_${chr}.fa.gz"

echo "Chromosome ${chr}"
echo "   HAPS     : ${haps_file}"
echo "   SAMPLE   : ${sample_file}"
echo "   MASK     : ${mask_file}"
echo "   ANCESTOR : ${ancestor_file}"

"${PREP_SCRIPT}" \
  --haps   "${haps_file}" \
  --sample "${sample_file}" \
  --ancestor "${ancestor_file}" \
  --mask     "${mask_file}" \
  -o "Archaic_African_masked_${chr}"

rc=$?
if [ $rc -eq 0 ]; then
  echo "Finished chromosome ${chr}"
else
  echo "Failed chromosome ${chr} (exit ${rc})"
fi

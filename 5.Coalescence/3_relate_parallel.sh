#!/bin/bash
#SBATCH --job-name=relate_parallel
#SBATCH --array=1-22                  
#SBATCH --cpus-per-task=2             
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=relate_parallel_chr%a.log

RELATE_PARALLEL="/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/scripts/RelateParallel/RelateParallel.sh"

input_prefix="Archaic_African_masked"
output_prefix="relate_Archaic_African"
coal_file="/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/African/Relate_input_files/coal_rates/1000G_auto.coal"
mu=4e-9         
N=30000         
memory=5        

chr=${SLURM_ARRAY_TASK_ID}
map_file="/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/African/Relate_input_files/GRCh38/recomb_map/genetic_map_chr${chr}.txt"

echo "Running RelateParallel for chr${chr}"
echo "   --haps   ${input_prefix}_${chr}.haps.gz"
echo "   --sample ${input_prefix}_${chr}.sample.gz"
echo "   --dist   ${input_prefix}_${chr}.dist.gz"
echo "   --map    ${map_file}"
echo "   --coal   ${coal_file}"
echo "   -m       ${mu}"
echo "   -N       ${N}"
echo "   --memory ${memory}"

"${RELATE_PARALLEL}" \
    --threads ${SLURM_CPUS_PER_TASK} \
    --haps    "${input_prefix}_${chr}.haps.gz" \
    --sample  "${input_prefix}_${chr}.sample.gz" \
    --dist    "${input_prefix}_${chr}.dist.gz" \
    --map     "${map_file}" \
    -N        ${N} \
    --coal    "${coal_file}" \
    -m        ${mu} \
    --memory  ${memory} \
    -o        "${output_prefix}_chr${chr}"

rc=$?
if [ $rc -ne 0 ]; then
    echo "RelateParallel failed on chr${chr} (exit ${rc})"
    exit $rc
fi

# gzip the two main outputs
gzip -f "${output_prefix}_chr${chr}.anc"
gzip -f "${output_prefix}_chr${chr}.mut"

echo "Chromosome ${chr} complete"

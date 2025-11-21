#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=4:00:00

# Loop over chromosomes 1 to 22 and convert phased VCFs to Relate format

for chr in {1..22}; do
    input_vcf="/lisc/scratch/admixlab/aigerim/shapeit5/Archaic/vcf/Archaic_African_chr${chr}"
    haps_out="haps/Archaic_African_chr${chr}.haps"
    sample_out="sample/Archaic_African_${chr}.sample"
    echo "Processing chromosome $chr"

    /lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/bin/RelateFileFormats --mode ConvertFromVcf \
        --haps "$haps_out" \
        --sample "$sample_out" \
        -i "$input_vcf"

    if [ $? -eq 0 ]; then
        echo "✅ Finished chromosome $chr"
    else
        echo "❌ Failed chromosome $chr"
    fi
done

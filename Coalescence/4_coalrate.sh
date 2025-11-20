#!/bin/bash
#SBATCH --job-name=coalrate
#SBATCH --mem=100G
#SBATCH --time=2:00:00

#nomask
/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/African/CoalRate/bin/CoalRate --mode local_ancestry --input relate_Archaic_African --output relate_AA_all_chr_noMaskHuman --chr all_chr.txt --poplabels all.inds.poplabels.space_sep --bins 2,7,0.2 --years_per_gen 30
#truemask
/lisc/scratch/admixlab/aigerim/shapeit5/African/vcf/relate/African/CoalRate/bin/CoalRate --mode local_ancestry --input relate_Archaic_African --output relate_AA_all_chr_MaskHuman --chr all_chr.txt --poplabels assignment_chrom_all.bed --bins 2,7,0.2 --years_per_gen 30

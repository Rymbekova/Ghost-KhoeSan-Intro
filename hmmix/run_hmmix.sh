#!/bin/bash
#SBATCH -N 2
#SBATCH --time=3:00:00
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --mem=100G

module load conda
conda activate sstar-analysis
export LD_LIBRARY_PATH=~/.conda/envs/sstar-analysis/lib/libpython3.8.so.1.0
export OPENBLAS_NUM_THREADS=64


while IFS= read -r foldername
do
    # Change the working directory to the current foldername
    cd "$foldername"

    # Loop through files foldername.1.ref50.tgt1.recode.bcf to foldername.23.ref50.tgt1.recode.bcf
    for i in {1..23}
    do
        sample="$foldername.$i.ref50.tgt1.recode"
        outgroup_file="outgroup.$foldername.$i.txt"
        mutationrate_file="mutationrate.$foldername.$i.bed"
        obs_file="obs.$i.$foldername.txt"

        # Perform hmmix operations for each sample
        hmmix create_outgroup -ind=individuals.json  -vcf="$sample.bcf" -out="$outgroup_file" 
        hmmix mutation_rate -outgroup="$outgroup_file" -window_size=1000000  -out "$mutationrate_file"
        hmmix create_ingroup -ind=individuals.json -vcf="$sample.bcf"  -out="obs.$i" -outgroup="$outgroup_file" 
        hmmix train -obs="$obs_file"  -mutrates="$mutationrate_file" -out=trained_"$sample".json
        hmmix decode -obs="$obs_file" -mutrates="$mutationrate_file" -param=trained_"$sample".json -out="$sample"
    done

    # Clean up files within the folder
    rm *.bed
    rm outgroup*
    rm obs*
    rm trained*

    # Move back to the parent directory
    cd ..

    echo "Processing for $foldername complete."
done < all.San.ind


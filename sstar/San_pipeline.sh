#!/bin/bash
#SBATCH --job-name=San_pipeline
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=12:00:00

# Check if the input file exists
input_file="all.San.ind"

if [ ! -f "$input_file" ]; then
    echo "Input file $input_file not found."
    exit 1
fi

# List of reference names
ref_names=(
    "HG02568"
    "HG02922"
    "HG03052"
    "HGDP_DNK02"
    "HGDP_HGDP00927"
    "HGDP_HGDP01284"
    "NA19017"
    "PNP010"
    "PNP011"
    "PNP012"
    "PNP013"
    "PNP014"
    "PNP030"
    "PNP031"
    "PNP032"
    "PNP033"
    "PNP034"
    "PNP050"
    "PNP051"
    "PNP052"
    "PNP053"
    "PNP054"
    "PNP080"
    "PNP081"
    "PNP082"
    "PNP083"
    "PNP084"
    "SGDP_LP6005441-DNA_E07"
    "SGDP_LP6005441-DNA_F07"
    "SGDP_LP6005442-DNA_A02"
    "SGDP_LP6005442-DNA_A10"
    "SGDP_LP6005442-DNA_B02"
    "SGDP_LP6005442-DNA_B10"
    "SGDP_LP6005442-DNA_D09"
    "SGDP_LP6005442-DNA_E11"
    "SGDP_LP6005442-DNA_F09"
    "SGDP_LP6005442-DNA_F11"
    "SGDP_LP6005442-DNA_G10"
    "SGDP_LP6005442-DNA_G11"
    "SGDP_LP6005442-DNA_H10"
    "SGDP_LP6005442-DNA_H11"
    "SGDP_LP6005443-DNA_B09"
    "SGDP_LP6005443-DNA_H08"
    "SGDP_SS6004480"
    "SGDP_SS6004475"
    "SGDP_SS6004470"
    "SGDP_LP6005677-DNA_G01"
    "SGDP_LP6005619-DNA_B01"
    "SGDP_LP6005619-DNA_C01"
    "SGDP_LP6005443-DNA_F06"
)

# Read sample names from the input file and create folders and files
while IFS= read -r sample_name; do
    # Create folder
    mkdir "$sample_name"
    
    # Create tgt.ind file
    echo "$sample_name" > "$sample_name/tgt.ind"
    
    # Create ref.ind file with the list of reference names
    for ref_name in "${ref_names[@]}"; do
        echo "$ref_name" >> "$sample_name/ref.ind"
    done
    
    # Create together.ind file by concatenating tgt.ind and ref.ind
    cat "$sample_name/tgt.ind" "$sample_name/ref.ind" > "$sample_name/together.ind"
    
    echo "Files for $sample_name created successfully."
    
    # Navigate to the folder
    cd "$sample_name"
    
    # Run the vcftools command with custom output name
    vcftools --vcf ../25KS.48RHG.74comp.HCBP.21.recalSNP99.9.recalINDEL99.0.vcf --keep together.ind --remove-indels --max-missing 1.0 --min-alleles 2 --max-alleles 2 --recode --out "$sample_name.ref50.tgt1"
    
    # Return to the parent directory
    cd ..
    
    echo "vcftools command for $sample_name executed successfully."
done < "$input_file"

echo "All folders and vcftools commands executed successfully."


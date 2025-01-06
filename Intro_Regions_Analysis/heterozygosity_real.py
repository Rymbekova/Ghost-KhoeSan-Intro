import sys
import subprocess

# Function to extract genotypes from a specific region in the VCF for a given individual
def extract_genotypes(region_vcf_individual):
    region, vcf_file, individual = region_vcf_individual
    command = f"bcftools view -m2 -M2 -v snps -r {region} {vcf_file} | bcftools query -s {individual} -f '%CHROM %POS %REF %ALT [ %GT ]\\n'"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return None
    return result.stdout.strip().split('\n')

# Function to calculate heterozygosity
def calculate_heterozygosity(genotypes, region_length):
    het_genotypes = sum(1 for gt in genotypes if gt in ['0/1', '1/0'])
    return (het_genotypes / region_length) * 1000 if region_length > 0 else 0

# Function to calculate the region length based on start and end positions
def calculate_region_length(region):
    chrom, positions = region.split(':')
    start, end = map(int, positions.split('-'))
    return end - start

# Main function to process the real region, VCF file, and individual
def process_real_region(region, vcf_file, individual):
    region_length = calculate_region_length(region)  # Calculate the region length
    genotypes_data = extract_genotypes((region, vcf_file, individual))

    if genotypes_data:
        genotypes = []
        for gt_line in genotypes_data:
            try:
                gt_columns = gt_line.split()[4:]  # Extract genotype columns
                genotypes.extend(gt_columns)      # Append all genotypes
            except ValueError:
                continue
        
        valid_genotypes = [gt for gt in genotypes if gt in ['0/0', '0/1', '1/0', '1/1']]
        
        if valid_genotypes:
            heterozygosity = calculate_heterozygosity(valid_genotypes, region_length)
            return heterozygosity  # Return the heterozygosity value
    return 0  # Return 0 if no valid genotypes or an error occurred

if __name__ == "__main__":
    # Capture the arguments passed to the script
    region = sys.argv[1]
    vcf_file = sys.argv[2]
    individual = sys.argv[3]  # Added individual argument
    
    # Process the real region for the specified individual
    heterozygosity_value = process_real_region(region, vcf_file, individual)
    
    # Print only the heterozygosity value
    print(heterozygosity_value)

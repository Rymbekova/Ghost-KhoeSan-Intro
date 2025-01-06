# Set custom library path
.libPaths("/lisc/scratch/admixlab/aigerim/Rlibs")

# Load necessary libraries
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)
library(parallel)

# Define paths and files
linsight_dir <- "/lisc/scratch/admixlab/aigerim/sstar-analysis-main/African/San21/linsight_scores/"
output_file <- "individual_San_mean_linsight_scores.csv"
vcf_template <- "/lisc/scratch/admixlab/aigerim/African/25KS.48RHG.74comp.HCBP.%d.recalSNP99.9.recalINDEL99.0.vcf.gz"

# Read individual names from 'all.San.ind' file
individuals <- readLines("all.San.ind")

# Initialize the output file with headers
write.csv(data.frame(Individual = character(), Mean_Real = numeric(), Mean_Random = numeric()),
          file = output_file, row.names = FALSE, quote = FALSE)

# Function to load each RDS file and convert it to GRanges with dequantized scores
load_linsight_as_granges <- function(chr) {
  rds_file <- file.path(linsight_dir, paste0("linsight.UCSC.hg19.", chr, ".rds"))
  if (file.exists(rds_file)) {
    linsight_chr <- readRDS(rds_file)
    starts <- cumsum(c(1, head(runLength(linsight_chr), -1)))
    ends <- cumsum(runLength(linsight_chr))
    scores <- (as.numeric(runValue(linsight_chr)) - 1) / 100
    GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends), score = scores)
  } else {
    stop("RDS file not found for chromosome: ", chr)
  }
}

# Calculate weighted mean LINSIGHT score for a specific region
calculate_weighted_mean <- function(region, linsight_data) {
  overlap_scores <- subsetByOverlaps(linsight_data, region)
  widths <- width(overlap_scores)
  scores <- mcols(overlap_scores)$score
  if (length(scores) > 0) {
    weighted_sum <- sum(scores * widths, na.rm = TRUE)
    total_width <- sum(widths, na.rm = TRUE)
    list(weighted_sum = weighted_sum, total_width = total_width)
  } else {
    list(weighted_sum = 0, total_width = 0)
  }
}

# Process each region independently in parallel
process_region <- function(region, individual, vcf_template) {
  tryCatch({
    # Identify chromosome and load corresponding LINSIGHT data
    chr <- as.character(seqnames(region))
    linsight_data <- load_linsight_as_granges(chr)
    
    # Define the region string for VCF extraction
    region_str <- paste0(chr, ":", start(region), "-", end(region))
    vcf_path <- sprintf(vcf_template, as.integer(sub("chr", "", chr)))
    
    # Use bcftools to extract specific region for individual from VCF
    temp_vcf <- tempfile(fileext = ".vcf.gz")
    system(paste("bcftools view -r", shQuote(region_str), "-s", shQuote(individual), "-Oz",
                 shQuote(vcf_path), "-o", shQuote(temp_vcf)))
    
    if (file.exists(temp_vcf) && file.info(temp_vcf)$size > 0) {
      vcf_data <- readVcf(temp_vcf, genome = "hg19")
      file.remove(temp_vcf)
      
      # Calculate weighted mean for the region
      mean_result <- calculate_weighted_mean(region, linsight_data)
      return(mean_result)
    } else {
      cat("No VCF data available for region", region_str, "\n")
      return(list(weighted_sum = 0, total_width = 0))
    }
  }, error = function(e) {
    cat("Error in processing region", region, "for individual", individual, ":", e$message, "\n")
    return(list(weighted_sum = 0, total_width = 0))
  })
}

# Process each individual
for (individual in individuals) {
  cat("Processing individual:", individual, "\n")
  
  # Define paths to BED files for real and random regions
  real_bed_path <- file.path(individual, paste0(individual, ".all.merged.tracts.bed"))
  random_bed_path <- file.path(individual, paste0(individual, ".random.merged.tracts.bed"))
  
  # Check if BED files exist; if not, skip this individual
  if (!file.exists(real_bed_path) || !file.exists(random_bed_path)) {
    cat("Skipping individual:", individual, "- BED files not found\n")
    next
  }
  
  # Import BED files as GRanges
  real_regions <- import(real_bed_path, format = "BED")
  random_regions <- import(random_bed_path, format = "BED")
  
  # Process real regions in parallel
  real_results <- mclapply(seq_along(real_regions), function(i) {
    process_region(real_regions[i], individual, vcf_template)
  }, mc.cores = 4)  # Adjust `mc.cores` as needed
  
  # Process random regions in parallel
  random_results <- mclapply(seq_along(random_regions), function(i) {
    process_region(random_regions[i], individual, vcf_template)
  }, mc.cores = 4)  # Adjust `mc.cores` as needed
  
  # Aggregate results
  cumulative_weighted_sum_real <- sum(sapply(real_results, function(x) x$weighted_sum))
  total_width_real <- sum(sapply(real_results, function(x) x$total_width))
  cumulative_weighted_sum_random <- sum(sapply(random_results, function(x) x$weighted_sum))
  total_width_random <- sum(sapply(random_results, function(x) x$total_width))
  
  # Calculate mean LINSIGHT scores for the individual
  mean_score_real <- if (total_width_real > 0) cumulative_weighted_sum_real / total_width_real else NA
  mean_score_random <- if (total_width_random > 0) cumulative_weighted_sum_random / total_width_random else NA
  
  # Append results to CSV file
  write.table(data.frame(Individual = individual, Mean_Real = mean_score_real, Mean_Random = mean_score_random),
              file = output_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  
  # Clean up memory for this individual
  rm(real_regions, random_regions, real_results, random_results)
  gc()
}

cat("Results saved to", output_file, "\n")

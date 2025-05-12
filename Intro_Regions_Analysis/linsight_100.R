# Set custom library path
.libPaths("/lisc/scratch/admixlab/aigerim/Rlibs")

# Load necessary libraries
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)
library(parallel)

# -----------------------------------------------------------------------------
# CONFIGURATION ----------------------------------------------------------------
# -----------------------------------------------------------------------------

linsight_dir   <- "/lisc/scratch/admixlab/aigerim/sstar-analysis-main/African/San21/linsight_scores/"
# New VCF template uses {chr} placeholder (e.g., 1, 2, X)
vcf_template   <- "/lisc/scratch/admixlab/aigerim/ensembl-vep/sstar_subset_filter1_%s.vcf.gz"

# A text file that lists all individuals (one per line)
individuals    <- readLines("all.San.ind")

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS -------------------------------------------------------------
# -----------------------------------------------------------------------------

# Load the LINSIGHT run‑length‑encoded RDS for a chromosome and convert to
# GRanges with de‑quantised scores.
load_linsight_as_granges <- function(chr) {
  rds_file <- file.path(linsight_dir, sprintf("linsight.UCSC.hg19.%s.rds", chr))
  if (!file.exists(rds_file)) {
    stop("LINSIGHT RDS file not found for chromosome ", chr)
  }
  rl <- readRDS(rds_file)
  starts <- cumsum(c(1L, head(runLength(rl), -1L)))
  ends   <- cumsum(runLength(rl))
  scores <- (as.numeric(runValue(rl)) - 1) / 100      # de‑quantise
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends), score = scores)
}

# Weighted mean LINSIGHT score for a single GRanges region
calculate_weighted_mean <- function(region, linsight_data) {
  hits   <- subsetByOverlaps(linsight_data, region)
  w      <- width(hits)
  s      <- mcols(hits)$score
  if (length(s) == 0) {
    return(list(sum = 0, width = 0))
  }
  list(sum = sum(s * w, na.rm = TRUE), width = sum(w, na.rm = TRUE))
}

# Process every interval within a BED for a single individual; returns the
# overall mean for that BED (real or a random replicate)
process_bed <- function(bed_path, individual, mc.cores = 4) {
  if (!file.exists(bed_path)) {
    warning("BED file not found: ", bed_path)
    return(NA_real_)
  }
  regions <- import(bed_path, format = "BED")
  if (length(regions) == 0) {
    warning("No regions in BED: ", bed_path)
    return(NA_real_)
  }

  # Split by chromosome so we only load each LINSIGHT RDS once per chromosome
  chrs <- unique(as.character(seqnames(regions)))
  linsight_cache <- setNames(vector("list", length(chrs)), chrs)

  means <- mclapply(seq_along(regions), function(i) {
    region <- regions[i]
    chr    <- as.character(seqnames(region))

    # Load & cache LINSIGHT for this chromosome
    if (is.null(linsight_cache[[chr]])) {
      linsight_cache[[chr]] <<- load_linsight_as_granges(chr)
    }
    lins_data <- linsight_cache[[chr]]

    # Extract overlapping VCF records for the individual
    region_str <- sprintf("%s:%d-%d", chr, start(region), end(region))
    # Remove "chr" prefix for VCF file names (adjust if your files include it)
    vcf_path   <- sprintf(vcf_template, sub("^chr", "", chr))
    tmp_vcf    <- tempfile(fileext = ".vcf.gz")

    system2("bcftools", c("view", "-r", region_str,
                           "-s", individual, "-Oz", vcf_path, "-o", tmp_vcf),
            stdout = FALSE, stderr = FALSE)

    on.exit(unlink(tmp_vcf), add = TRUE)

    if (!file.exists(tmp_vcf) || file.info(tmp_vcf)$size == 0) {
      return(list(sum = 0, width = 0))
    }

    vcf <- readVcf(tmp_vcf, genome = "hg19")  # retained for completeness

    calculate_weighted_mean(region, lins_data)
  }, mc.cores = mc.cores)

  total_sum   <- sum(vapply(means, `[[`, numeric(1), "sum"))
  total_width <- sum(vapply(means, `[[`, numeric(1), "width"))
  if (total_width == 0) NA_real_ else total_sum / total_width
}

# -----------------------------------------------------------------------------
# MAIN LOOP --------------------------------------------------------------------
# -----------------------------------------------------------------------------

for (ind in individuals) {
  message("\n>>> Processing ", ind)

  # Real BED
  real_bed <- file.path(ind, sprintf("%s_overlap.bed", ind))
  real_mean <- process_bed(real_bed, ind)

  # 100 random replicates
  random_means <- numeric(100)
  for (rep in 1:100) {
    random_bed <- file.path(ind, sprintf("%s.random.replicate%d_overlap.bed", ind, rep))
    random_means[rep] <- process_bed(random_bed, ind)
  }

  # Build output data frame (one row, 102 columns: Individual + Real + 100 Randoms)
  out <- data.frame(Individual = ind,
                    Real       = real_mean,
                    setNames(as.data.frame(t(random_means)),
                             sprintf("Random_%d", 1:100)))

  # Write individual‑specific CSV
  out_file <- sprintf("%s_linsight_100.csv", ind)
  write.csv(out, file = out_file, row.names = FALSE, quote = FALSE)
  message("    → Saved results to ", out_file)

  # Tidy‑up cached LINSIGHT objects to free memory between individuals
  gc()
}

message("\nAll done!")

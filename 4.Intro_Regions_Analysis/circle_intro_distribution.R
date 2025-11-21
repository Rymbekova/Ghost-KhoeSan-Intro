# Install necessary packages
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

# Load the packages
library(rtracklayer)
library(circlize)
library(data.table)

# Define the path to your BED files
bed_files <- list.files(path = "/home/rymbekovaa95/Desktop/African/", pattern = "*overlap.bed", full.names = TRUE)

# Function to read a BED file and return a data frame
read_bed_file <- function(file) {
  bed <- import(file)
  as.data.frame(bed)
}

# Read all BED files and store in a list
bed_list <- lapply(bed_files, read_bed_file)

# Check the structure of one of the BED files to ensure it's read correctly
print(head(bed_list[[1]]))

# Function to aggregate data in sliding windows
aggregate_in_sliding_windows <- function(bed_list, window_size = 2e6, step_size = 2e6) {
  combined <- do.call(rbind, bed_list)
  combined <- combined[order(combined$seqnames, combined$start), ]
  
  # Create an empty data table for aggregated results
  aggregated <- data.table(chr = character(), start = numeric(), end = numeric(), value = numeric())
  
  # Process each chromosome separately
  chromosomes <- unique(combined$seqnames)
  for (chr in chromosomes) {
    chr_data <- combined[combined$seqnames == chr, ]
    max_pos <- max(chr_data$end)
    for (start_pos in seq(0, max_pos, by = step_size)) {
      end_pos <- start_pos + window_size
      in_window <- chr_data[chr_data$start < end_pos & chr_data$end > start_pos, ]
      total_value <- sum(in_window$width)
      aggregated <- rbind(aggregated, data.table(chr = chr, start = start_pos, end = end_pos, value = total_value))
    }
  }
  
  aggregated
}

# Aggregate the data
aggregated_data <- aggregate_in_sliding_windows(bed_list)

# Check the structure of the aggregated data
print(head(aggregated_data))

# Save the plot as PDF
#pdf("line_circular_plot.pdf", width = 8, height = 8)

# Initialize the circular plot
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))

# Add aggregated data as inner blue circle
circos.genomicTrackPlotRegion(aggregated_data, ylim = c(0, max(aggregated_data$value)),
                              track.height = 0.3, panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, col = "blue", ...)
                              })

# Clear the plot
#circos.clear()

# Close the PDF device
#dev.off()
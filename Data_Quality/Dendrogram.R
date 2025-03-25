personal_lib <- "~/R/library"

if (!dir.exists(personal_lib)) {
    dir.create(personal_lib, recursive = TRUE)
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

library(ggplot2)

vcf.fn <- "25KS.48RHG.74comp.HCBP.21.recalSNP99.9.recalINDEL99.0.vcf.gz"
vcf_data <- read.vcf(vcf.fn)  

genotypes <- as.data.frame(vcf_data)
genotype_matrix <- as.matrix(genotypes)

distance_matrix <- dist(genotype_matrix)
hc <- hclust(distance_matrix, method = "average")

pdf("dendrogram_plot_all.pdf")
plot(hc, main = "Dendrogram based on Genetic Distance")
dev.off()

ggsave(filename = "dendrogram_plot_all.pdf", plot = last_plot())

savehistory(file = "dendrogram.R")


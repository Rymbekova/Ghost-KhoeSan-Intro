#!/usr/bin/env Rscript

library(ggplot2)
require(gridExtra)

options(scipen=100)
'%ni%' <- Negate('%in%')

mura <- 1.29e-08  # mutation rate
gtim <- 30        # generation time

# Function to calculate standard deviation across matrices
sdfun <- function(inputtab) {
  vec <- unlist(inputtab, use.names = FALSE)
  DIM <- dim(inputtab[[1]])
  n <- length(inputtab)
  list.sd <- tapply(vec, rep(1:prod(DIM), times = n), sd)
  attr(list.sd, "dim") <- DIM
  as.data.frame(list.sd)
}

# Function to compute admixture time and coalescence
tfun <- function(tb) {
  p <- tb[3, 2]
  q <- tb[4, 1]
  Tadm <- (p + (2 * q)) / (2 * 1000 * mura)
  c(Tadm * gtim, (p / (Tadm * 2 * 1000 * mura)))
}

# Wrapper to extract useful values
colfun <- function(input) {
  c(tfun(input), 
    round(input[2, 2] / (1000 * mura) * gtim), 
    round(input[2, 1] / (1000 * mura) * gtim))
}

# Load individual sample list
inds <- list()
inds[[1]] <- unlist(read.table("all.San.ind", as.is = TRUE))
gn <- "San"

# Read .hmm files and extract matrices
tabls <- list()
itab <- list()
tabls[[1]] <- list()
itab[[1]] <- list()

for (i in 1:length(inds[[1]])) {
  hmmfile <- paste0("averaged_", inds[[1]][i], ".hmm")
  if (!file.exists(hmmfile)) next
  t1 <- try(readLines(hmmfile), silent=TRUE)
  if (inherits(t1, "try-error")) next
  
  # Extract parameters safely using regex
  start_line <- grep("^starting", t1, value = TRUE)
  emis_line <- grep("^emissions", t1, value = TRUE)
  trans_line <- grep("^transitions", t1, value = TRUE)
  
  start_probs <- as.numeric(unlist(regmatches(start_line, gregexpr("[0-9.]+", start_line))))
  emissions   <- as.numeric(unlist(regmatches(emis_line, gregexpr("[0-9.]+", emis_line))))
  trans_vals  <- as.numeric(unlist(regmatches(trans_line, gregexpr("[0-9.]+", trans_line))))
  
  # Check length before continuing
  if (length(start_probs) != 2 || length(emissions) != 2 || length(trans_vals) != 4) {
    warning(paste("Skipping", hmmfile, "- invalid parameter lengths"))
    next
  }
  
  transitions <- matrix(trans_vals, ncol = 2, byrow = TRUE)
  mat <- rbind(start_probs, emissions, transitions)
  tabls[[1]][[i]] <- mat
  itab[[1]][[i]] <- t1
}

# Summary stats
reman <- list()
calval <- list()
calval2 <- list()
crval <- list()
crval2 <- list()
coval <- list()

t2 <- Reduce("+", tabls[[1]]) / length(tabls[[1]])
rownames(t2) <- c("Starting", "Emission", "Transition_within", "Transition_arch")
reman[[1]] <- t2

sdv <- sdfun(tabls[[1]]) * 2.58
crval[[1]] <- do.call(rbind, lapply(tabls[[1]], tfun))

crval2[[1]] <- c(colMeans(crval[[1]]) * c(1, 100), 
                 sd(crval[[1]][, 1]) * 2.58, 
                 sd(crval[[1]][, 2]) * 2.58 * 100)

calval[[1]] <- c(round(t2[1,2]*100,1), 
                 round(t2[2,2] / (1000 * mura) * gtim),
                 round(t2[2,1] / (1000 * mura) * gtim))

calval2[[1]] <- c(round((t2[1,2] - sdv[1,2]) * 100, 1),
                  round((t2[1,2] + sdv[1,2]) * 100, 1),
                  round((t2[2,2] - sdv[2,2]) / (1000 * mura) * gtim),
                  round((t2[2,2] + sdv[2,2]) / (1000 * mura) * gtim),
                  round((t2[2,1] - sdv[2,1]) / (1000 * mura) * gtim),
                  round((t2[2,1] + sdv[2,1]) / (1000 * mura) * gtim))

coval[[1]] <- do.call(rbind, lapply(tabls[[1]], colfun))

crrval <- cbind(inds[[1]], coval[[1]])
crrval
colnames(crrval) <- c("", "T_admix", "P_admix", "Tcoal_extern", "Tcoal_intern")
write.table(crrval, file=paste0("mut_gen_", gn, ".tsv"), sep="\t", row.names=FALSE, quote=FALSE)

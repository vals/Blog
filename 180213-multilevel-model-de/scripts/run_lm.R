#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop('Exactly 3 arguments required (metadata, genedata, outfile)')
}

metadata_path <- args[1]
genedata_path <- args[2]
outfile_path <- args[3]

library(tidyverse)
library(lmtest)

metadata <- read_csv(metadata_path)
gene_data <- read_csv(genedata_path)

# Define models
f1 <- 'log(expr + 1) ~ log(total_count)'
f2 <- 'log(expr + 1) ~ log(total_count) + stim'

# Function for fitting the two models to each gene and doing an LR test.

fit_lm <- function(expr_vec) {
  tmp_df <- cbind(metadata, data.frame(expr = unname(expr_vec)))
  
  mres1 <- lm(f1, data = tmp_df)
  mres2 <- lm(f2, data = tmp_df)
  
  lrt_res <- lrtest(mres1, mres2)
  pval <- lrt_res$`Pr(>Chisq)`[2]
  effect_size <- unname(coef(mres2)[3])
  
  c(pval = pval, effect_size = effect_size)
}

e_fit_lm <- function(expr_vec) {
  tryCatch(fit_lm(expr_vec), error = function(e) c(pval = NA, effect_size = NA))
}

# Apply the function

fit_results <- apply(gene_data, 2, e_fit_lm)

# Format results

results <- as.tibble(cbind(t(fit_results), data.frame(gene = rownames(t(fit_results)))))

write_csv(results, outfile_path)

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop('Exactly 3 arguments required (metadata, genedata, outfile)')
}

metadata_path <- args[1]
genedata_path <- args[2]
outfile_path <- args[3]

library(tidyverse)
library(lme4)

metadata <- read_csv(metadata_path)
gene_data <- read_csv(genedata_path)

# Define models
f1 <- 'expr ~ offset(log(total_count)) + (1 | sample)'
f2 <- 'expr ~ offset(log(total_count)) + (1 | stim/sample)'


# Function for fitting the two models to each gene and doing an LR test.

fit_lme <- function(expr_vec) {
  tmp_df <- cbind(metadata, data.frame(expr = unname(expr_vec)))
  
  mres1 <- glmer(f1, data = tmp_df, family = poisson())
  mres2 <- glmer(f2, data = tmp_df, family = poisson())
  anova_res <- anova(mres1, mres2)
  pval <- anova_res$`Pr(>Chisq)`[2]
  effect_size <- ranef(mres2)$stim[,1][1] - ranef(mres2)$stim[,1][2]
  
  c(pval = pval, effect_size = effect_size)
}

e_fit_lme <- function(expr_vec) {
  tryCatch(fit_lme(expr_vec), error = function(e) c(pval = NA, effect_size = NA))
}

# Apply the function

fit_results <- apply(gene_data, 2, e_fit_lme)

# Format results

results <- as.tibble(cbind(t(fit_results), data.frame(gene = rownames(t(fit_results)))))

write_csv(results, outfile_path)


# qplot(results$effect_size, -log10(results$pval))

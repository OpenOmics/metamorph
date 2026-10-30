#!/usr/bin/env Rscript
# Called by expand_taxonomy.py; do not feed this worker overlapping parent rows.
suppressPackageStartupMessages({
  library(phyloseq)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, '--check')) quit(status = 0)
if (length(args) != 2L) stop('Expected species.csv output_directory')
ranks <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
x <- fread(args[1], check.names = FALSE)
if (!identical(names(x)[seq_len(7)], ranks)) stop('Invalid taxonomy columns')
samples <- names(x)[-(seq_len(7))]
tax <- as.matrix(x[, ..ranks])
abundance <- as.matrix(x[, ..samples])
storage.mode(abundance) <- 'double'
ids <- paste0('feature_', seq_len(nrow(x)))
rownames(tax) <- rownames(abundance) <- ids
ps <- phyloseq(otu_table(abundance, taxa_are_rows = TRUE), tax_table(tax))
expected <- sample_sums(ps)
rank_totals <- list()
# Aggregate progressively from species toward kingdom. Every intermediate rank
# has a nonempty placeholder when unresolved, so no abundance is dropped.
for (i in rev(seq_along(ranks))) {
  if (i < length(ranks)) ps <- tax_glom(ps, taxrank = ranks[i], NArm = FALSE)
  mat <- as(otu_table(ps), 'matrix')
  if (!taxa_are_rows(ps)) mat <- t(mat)
  tx <- as(tax_table(ps), 'matrix')
  paths <- apply(tx[, seq_len(i), drop = FALSE], 1, paste, collapse = '|')
  if (anyDuplicated(paths)) stop('Duplicate lineages after aggregation')
  if (any(abs(colSums(mat) - expected) > 1e-7 * pmax(1, abs(expected)))) {
    stop('Abundance conservation failed at ', ranks[i])
  }
  tab <- as.data.table(mat)
  setnames(tab, colnames(mat))
  tab[, clade_name := paths]
  setcolorder(tab, c('clade_name', samples))
  setorder(tab, clade_name)
  fwrite(tab, file.path(args[2], paste0(tolower(ranks[i]), '.csv')), quote = 'auto')
  fwrite(tab, file.path(args[2], 'all_ranks.csv'), quote = 'auto', append = i < length(ranks))
  rank_totals[[i]] <- data.table(sample = names(expected), rank = ranks[i],
                                total = as.numeric(colSums(mat)),
                                within_100_tolerance = abs(colSums(mat) - 100) <= 1e-7)
}
fwrite(rbindlist(rank_totals), file.path(args[2], 'rank_totals.csv'))
writeLines(capture.output(sessionInfo()), file.path(args[2], 'R_session.txt'))

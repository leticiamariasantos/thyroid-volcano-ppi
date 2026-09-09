# Multi-compartment contribution to variance of top raw DEGs.
suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })
source(here("scripts", "phase2", "20_upgrade_config.R"))

expr <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"))
meta <- fread(file.path(DIR_VARIANTS, "analysis_metadata.tsv"))
comp <- readRDS(file.path(DIR_DECONV, "composition_consensus.rds"))
de <- fread(file.path(DIR_DE_MULTI, "raw", "voom_qw_full_results.tsv"))
genes <- head(de[order(adj.P.Val, -abs(logFC)), gene_symbol], 100L)
genes <- intersect(genes, rownames(expr))

comp <- as.data.frame(comp[match(meta$sample, rownames(comp)), , drop = FALSE])
names(comp) <- make.names(names(comp))
dat0 <- cbind(as.data.frame(meta), comp)
base_formula <- ~ condition + source
full_formula <- as.formula(paste("~ condition + source +", paste(names(comp), collapse = " + ")))

one_gene <- function(gene) {
  dat <- dat0; dat$y <- as.numeric(expr[gene, ])
  base <- lm(update(base_formula, y ~ .), dat)
  full <- lm(update(full_formula, y ~ .), dat)
  sse_base <- deviance(base); sse_full <- deviance(full)
  total <- pmax(0, 100 * (sse_base - sse_full) / sse_base)
  parts <- rbindlist(lapply(names(comp), function(cell) {
    reduced_terms <- setdiff(names(comp), cell)
    reduced_formula <- as.formula(paste("y ~ condition + source",
      if (length(reduced_terms)) paste("+", paste(reduced_terms, collapse = "+")) else ""))
    reduced <- lm(reduced_formula, dat)
    data.table(gene = gene, cell_type = cell,
      percent_variance_unique = pmax(0, pmin(100, 100 * (deviance(reduced) - sse_full) / sse_base)),
      percent_variance_all_composition = total)
  }))
  parts
}
contrib <- rbindlist(lapply(genes, one_gene))
contrib <- merge(contrib, de[, .(gene = gene_symbol, raw_logFC = logFC, raw_FDR = adj.P.Val)], by = "gene")
fwrite_tsv(contrib, file.path(DIR_DECONV, "top_deg_celltype_variance.tsv"))

top_plot <- contrib[, .(total = max(percent_variance_all_composition)), by = gene][order(-total)][1:min(30L, .N), gene]
p <- ggplot(contrib[gene %in% top_plot], aes(x = reorder(gene, percent_variance_all_composition),
                                             y = percent_variance_unique, fill = cell_type)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Variância única explicada (% do SSE base)", fill = "Compartimento",
       title = "Contribuição composicional nos top DEGs",
       subtitle = "R² parcial; condição e source permanecem no modelo") +
  theme_minimal(base_size = 10)
ggsave(file.path(DIR_FIG, "Fig23_Composition_Variance.png"), p, width = 9, height = 8, dpi = 300)
log_msg("Composition variance complete for", length(genes), "genes")

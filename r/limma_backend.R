#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: limma_backend.R expression.csv metadata.csv formula reference test output.csv")
}

suppressPackageStartupMessages(library(limma))

expression_path <- args[[1]]
metadata_path <- args[[2]]
design_formula <- args[[3]]
reference_level <- args[[4]]
test_level <- args[[5]]
output_path <- args[[6]]

expression <- read.csv(expression_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)
if (!identical(colnames(expression), rownames(metadata))) {
  if (!all(colnames(expression) %in% rownames(metadata))) {
    stop("Expression samples and metadata rows do not match.")
  }
  metadata <- metadata[colnames(expression), , drop = FALSE]
}

metadata$condition <- factor(metadata$condition, levels = c(reference_level, test_level))
for (column in setdiff(colnames(metadata), "condition")) {
  if (is.character(metadata[[column]])) metadata[[column]] <- factor(metadata[[column]])
}

design <- model.matrix(as.formula(design_formula), data = metadata)
if (qr(design)$rank < ncol(design)) stop("The limma design matrix is not full rank.")
condition_columns <- grep("^condition", colnames(design), value = TRUE)
if (length(condition_columns) != 1) {
  stop("Could not identify a unique test-versus-reference condition coefficient.")
}

fit <- lmFit(as.matrix(expression), design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
table <- topTable(fit, coef = condition_columns[[1]], number = Inf, sort.by = "none")

result <- data.frame(
  Probe_ID = rownames(table),
  MeanExpression = table$AveExpr,
  Log2FC = table$logFC,
  SE = abs(table$logFC / table$t),
  PValue = table$P.Value,
  FDR = table$adj.P.Val,
  Moderated_t = table$t,
  B_statistic = table$B,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(result, output_path, row.names = FALSE, na = "")

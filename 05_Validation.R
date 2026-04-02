# ============================================
# Step 5: Complete Runnable Validation Script
# ============================================
#
# Purpose:
#   1. Load prepared GSE75010 data and ML model (4 best genes: PROCR, FLT1, SPAG4, COL17A1)
#   2. Load and preprocess GSE25906 and GSE60438 raw expression data
#   3. Match genes by gene name (no GPL annotation required)
#   4. Validate ML model on independent datasets
#   5. Save validation results as RDS and CSV files
#
# Dependencies: GEOquery, limma, preprocessCore, pROC

# ============================================
# Configuration
# ============================================
FALLBACK_MAX_GENES       <- 20   # max candidates to scan in ML importance fallback
CLASSIFICATION_THRESHOLD <- 0.5  # probability cut-off for sensitivity/specificity

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║          Step 5: Validation on Independent Datasets          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

library(GEOquery)
library(limma)
library(preprocessCore)
library(pROC)

# ============================================
# [5.1] Load ML results
# ============================================

cat("[5.1] Loading ML results...\n\n")

if (!file.exists("results/04_ML_results.rds")) {
  stop("results/04_ML_results.rds not found. Please run Step 4 first.")
}

ml_results <- readRDS("results/04_ML_results.rds")
best_genes  <- ml_results$best_genes

cat("Best genes from ML model:\n")
for (i in seq_along(best_genes)) {
  cat(sprintf("  %d. %s\n", i, best_genes[i]))
}
cat(sprintf("\nTotal: %d genes\n\n", length(best_genes)))

# ============================================
# Helper: log2 normalise + quantile normalise
# ============================================

preprocess_expr <- function(expr_raw) {
  expr_log <- log2(expr_raw + 1)
  min_val   <- min(expr_log, na.rm = TRUE)
  if (min_val < 0) {
    expr_log <- expr_log - min_val + 1
  }
  # Remove rows where every value is missing (NA or NaN)
  nan_rows  <- rowSums(is.na(expr_log) | is.nan(expr_log)) == ncol(expr_log)
  expr_log  <- expr_log[!nan_rows, , drop = FALSE]
  expr_norm <- normalize.quantiles(as.matrix(expr_log))
  rownames(expr_norm) <- rownames(expr_log)
  colnames(expr_norm) <- colnames(expr_log)
  return(expr_norm)
}

# Helper: extract gene symbols (fall back to probe ID when Symbol is empty)
get_gene_symbols <- function(feature_data) {
  if ("Symbol" %in% colnames(feature_data)) {
    symbols <- as.character(feature_data$Symbol)
  } else {
    symbols <- rep("", nrow(feature_data))
  }
  missing  <- is.na(symbols) | symbols == ""
  if ("ID" %in% colnames(feature_data)) {
    symbols[missing] <- as.character(feature_data$ID[missing])
  } else {
    symbols[missing] <- as.character(seq_len(sum(missing)))
  }
  return(symbols)
}

# ============================================
# [5.2] Load GSE75010 (Discovery set) — reference gene names
# ============================================

cat("[5.2] Loading GSE75010 (Discovery set)...\n\n")

file_75010 <- "GSE75010_series_matrix.txt.gz"
if (!file.exists(file_75010)) {
  stop(paste("File not found:", file_75010))
}

eset_75010        <- getGEO(filename = file_75010, getGPL = FALSE)
expr_75010_raw    <- exprs(eset_75010)
feature_75010     <- fData(eset_75010)
gene_symbols_75010 <- get_gene_symbols(feature_75010)

cat("GSE75010 dimensions:", dim(expr_75010_raw), "\n")
cat("Best genes present in GSE75010:\n")
for (gene in best_genes) {
  status <- ifelse(gene %in% gene_symbols_75010, "✓", "✗")
  cat(sprintf("  %s %s\n", status, gene))
}
cat("\n")

# ============================================
# [5.3] Load GSE25906 (Validation set 1)
# ============================================

cat("[5.3] Loading GSE25906...\n\n")

file_25906 <- "GSE25906_series_matrix.txt.gz"
if (!file.exists(file_25906)) {
  stop(paste("File not found:", file_25906))
}

eset_25906         <- getGEO(filename = file_25906, getGPL = FALSE)
expr_25906_raw     <- exprs(eset_25906)
pheno_25906        <- pData(eset_25906)
feature_25906      <- fData(eset_25906)

expr_25906_norm    <- preprocess_expr(expr_25906_raw)
gene_symbols_25906 <- get_gene_symbols(feature_25906)

# Sample grouping — expects a "classification:ch1" column
if ("classification:ch1" %in% colnames(pheno_25906)) {
  classification_25906 <- pheno_25906$`classification:ch1`
} else {
  # Attempt fallback: look for title column
  classification_25906 <- pheno_25906$title
  cat("  ⚠ 'classification:ch1' not found; using 'title' for grouping.\n")
}
group_25906 <- factor(ifelse(grepl("preeclamp", classification_25906, ignore.case = TRUE),
                             "PE", "Control"),
                      levels = c("Control", "PE"))

cat("GSE25906 dimensions:", dim(expr_25906_norm), "\n")
cat("Samples - PE:", sum(group_25906 == "PE"),
    "  Control:", sum(group_25906 == "Control"), "\n")
cat("Best genes present in GSE25906:\n")
for (gene in best_genes) {
  status <- ifelse(gene %in% gene_symbols_25906, "✓", "✗")
  cat(sprintf("  %s %s\n", status, gene))
}
cat("\n")

# ============================================
# [5.4] Load GSE60438 (Validation set 2)
# ============================================

cat("[5.4] Loading GSE60438...\n\n")

# Accept either the combined or the GPL10558-specific matrix file
file_60438_candidates <- c("GSE60438-GPL10558_series_matrix.txt.gz",
                           "GSE60438_series_matrix.txt.gz")
file_60438 <- NULL
for (f in file_60438_candidates) {
  if (file.exists(f)) {
    file_60438 <- f
    break
  }
}
if (is.null(file_60438)) {
  stop(paste("GSE60438 series matrix not found. Expected one of:",
             paste(file_60438_candidates, collapse = ", ")))
}

eset_60438         <- getGEO(filename = file_60438, getGPL = FALSE)
expr_60438_raw     <- exprs(eset_60438)
pheno_60438        <- pData(eset_60438)
feature_60438      <- fData(eset_60438)

expr_60438_norm    <- preprocess_expr(expr_60438_raw)
gene_symbols_60438 <- get_gene_symbols(feature_60438)

# Sample grouping based on title (e.g. "PE_xxx" vs "Control_xxx")
title_60438 <- pheno_60438$title
group_60438 <- factor(ifelse(grepl("^PE", title_60438, ignore.case = TRUE),
                             "PE", "Control"),
                      levels = c("Control", "PE"))

cat("GSE60438 dimensions:", dim(expr_60438_norm), "\n")
cat("Samples - PE:", sum(group_60438 == "PE"),
    "  Control:", sum(group_60438 == "Control"), "\n")
cat("Best genes present in GSE60438:\n")
for (gene in best_genes) {
  status <- ifelse(gene %in% gene_symbols_60438, "✓", "✗")
  cat(sprintf("  %s %s\n", status, gene))
}
cat("\n")

# ============================================
# [5.5] Determine genes available in all datasets
# ============================================

cat("[5.5] Resolving gene availability...\n\n")

missing_25906 <- setdiff(best_genes, gene_symbols_25906)
missing_60438 <- setdiff(best_genes, gene_symbols_60438)

if (length(missing_25906) > 0) {
  cat("Genes missing in GSE25906:", paste(missing_25906, collapse = ", "), "\n")
}
if (length(missing_60438) > 0) {
  cat("Genes missing in GSE60438:", paste(missing_60438, collapse = ", "), "\n")
}

available_genes <- intersect(best_genes,
                             intersect(gene_symbols_25906, gene_symbols_60438))

cat(sprintf("Genes available in both validation sets: %d / %d\n",
            length(available_genes), length(best_genes)))
cat("Genes:", paste(available_genes, collapse = ", "), "\n\n")

# If fewer than 2 genes are shared, fall back to ML importance ranking
if (length(available_genes) < 2) {
  cat("⚠ Too few genes — falling back to ML importance ranking...\n\n")

  rf_importance   <- ml_results$rf_importance
  available_genes <- character(0)

  for (i in seq_len(min(FALLBACK_MAX_GENES, nrow(rf_importance)))) {
    gene <- rf_importance$Gene[i]
    if (gene %in% gene_symbols_25906 && gene %in% gene_symbols_60438) {
      available_genes <- c(available_genes, gene)
    }
  }

  cat(sprintf("Fallback: %d genes selected from ML importance ranking:\n",
              length(available_genes)))
  print(available_genes)
  cat("\n")

  if (length(available_genes) == 0) {
    stop("No genes found in either validation dataset. Cannot proceed.")
  }
}

# ============================================
# [5.6] Extract expression matrices (gene × sample)
# ============================================

cat("[5.6] Extracting expression matrices...\n\n")

extract_gene_expr <- function(expr_norm, gene_symbols, genes, dataset_name) {
  result_list <- list()
  for (gene in genes) {
    idx <- which(gene_symbols == gene)
    if (length(idx) > 0) {
      if (length(idx) > 1) {
        cat(sprintf("  ⚠ %s has %d probes in %s; using the first one.\n",
                    gene, length(idx), dataset_name))
      }
      result_list[[gene]] <- expr_norm[idx[1], ]
      cat(sprintf("  ✓ %s found in %s\n", gene, dataset_name))
    } else {
      cat(sprintf("  ✗ %s NOT found in %s\n", gene, dataset_name))
    }
  }
  if (length(result_list) == 0) {
    stop(paste("No genes extracted from", dataset_name))
  }
  mat <- do.call(rbind, result_list)
  rownames(mat) <- names(result_list)
  return(mat)
}

cat("GSE25906:\n")
expr_val_25906 <- extract_gene_expr(expr_25906_norm, gene_symbols_25906,
                                    available_genes, "GSE25906")
cat(sprintf("  → Matrix: %d genes × %d samples\n\n", nrow(expr_val_25906), ncol(expr_val_25906)))

cat("GSE60438:\n")
expr_val_60438 <- extract_gene_expr(expr_60438_norm, gene_symbols_60438,
                                    available_genes, "GSE60438")
cat(sprintf("  → Matrix: %d genes × %d samples\n\n", nrow(expr_val_60438), ncol(expr_val_60438)))

# ============================================
# [5.7] Validate: logistic regression + ROC
# ============================================

cat("[5.7] Building validation models...\n\n")

run_validation <- function(expr_data, group, dataset_name) {

  # Transpose to samples × genes; apply make.names() upfront for consistent names
  gene_cols <- make.names(rownames(expr_data))
  expr_mat  <- t(expr_data)
  data_val  <- as.data.frame(expr_mat)
  colnames(data_val) <- gene_cols
  data_val$PE_status <- as.numeric(group) - 1   # 0 = Control, 1 = PE
  formula_str  <- paste("PE_status ~", paste(gene_cols, collapse = " + "))

  # Guard against degenerate group composition
  n_positive <- sum(data_val$PE_status == 1)
  n_negative <- sum(data_val$PE_status == 0)
  if (n_positive == 0 || n_negative == 0) {
    cat(sprintf("✗ %s — Error: only one class present (PE=%d, Control=%d)\n\n",
                dataset_name, n_positive, n_negative))
    return(list(auc = NA, sensitivity = NA, specificity = NA,
                accuracy = NA, roc_obj = NULL, model = NULL, predictions = NULL))
  }

  tryCatch({
    model      <- glm(as.formula(formula_str), family = "binomial", data = data_val)
    pred_prob  <- predict(model, type = "response")
    roc_obj    <- roc(data_val$PE_status, pred_prob, quiet = TRUE)
    auc_val    <- as.numeric(roc_obj$auc)

    pred_class  <- ifelse(pred_prob > CLASSIFICATION_THRESHOLD, 1, 0)
    tp <- sum(pred_class == 1 & data_val$PE_status == 1)
    tn <- sum(pred_class == 0 & data_val$PE_status == 0)
    fp <- sum(pred_class == 1 & data_val$PE_status == 0)
    fn <- sum(pred_class == 0 & data_val$PE_status == 1)

    sensitivity <- tp / (tp + fn)   # n_positive > 0 guaranteed above
    specificity <- tn / (tn + fp)   # n_negative > 0 guaranteed above
    accuracy    <- (tp + tn) / (tp + tn + fp + fn)

    cat(sprintf("✓ %s\n", dataset_name))
    cat(sprintf("  - N samples  : %d\n", nrow(data_val)))
    cat(sprintf("  - AUC        : %.4f\n", auc_val))
    cat(sprintf("  - Sensitivity: %.3f\n", sensitivity))
    cat(sprintf("  - Specificity: %.3f\n", specificity))
    cat(sprintf("  - Accuracy   : %.3f\n\n", accuracy))

    return(list(
      auc         = auc_val,
      sensitivity = sensitivity,
      specificity = specificity,
      accuracy    = accuracy,
      roc_obj     = roc_obj,
      model       = model,
      predictions = pred_prob
    ))

  }, error = function(e) {
    cat(sprintf("✗ %s — Error: %s\n\n", dataset_name, e$message))
    return(list(auc = NA, sensitivity = NA, specificity = NA,
                accuracy = NA, roc_obj = NULL, model = NULL, predictions = NULL))
  })
}

result_25906 <- run_validation(expr_val_25906, group_25906, "GSE25906")
result_60438 <- run_validation(expr_val_60438, group_60438, "GSE60438")

# ============================================
# [5.8] Summary table
# ============================================

cat("[5.8] Validation summary...\n\n")

validation_results <- data.frame(
  Dataset    = c("GSE25906", "GSE60438"),
  N_Samples  = c(length(group_25906), length(group_60438)),
  N_PE       = c(sum(group_25906 == "PE"),      sum(group_60438 == "PE")),
  N_Control  = c(sum(group_25906 == "Control"), sum(group_60438 == "Control")),
  N_Genes    = c(nrow(expr_val_25906),          nrow(expr_val_60438)),
  AUC        = c(round(result_25906$auc, 4),    round(result_60438$auc, 4)),
  Sensitivity= c(round(result_25906$sensitivity, 3), round(result_60438$sensitivity, 3)),
  Specificity= c(round(result_25906$specificity, 3), round(result_60438$specificity, 3)),
  Accuracy   = c(round(result_25906$accuracy, 3),    round(result_60438$accuracy, 3)),
  stringsAsFactors = FALSE
)

cat("VALIDATION SUMMARY:\n")
cat(paste(rep("═", 62), collapse = ""), "\n\n")
print(validation_results)
cat("\n")

# ============================================
# [5.9] Save results
# ============================================

cat("[5.9] Saving results...\n\n")

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

validation_final <- list(
  best_genes         = best_genes,
  available_genes    = available_genes,
  validation_results = validation_results,
  result_25906       = result_25906,
  result_60438       = result_60438,
  expr_val_25906     = expr_val_25906,
  expr_val_60438     = expr_val_60438,
  group_25906        = group_25906,
  group_60438        = group_60438
)

saveRDS(validation_final, "results/05_Validation_Final.rds")
write.csv(validation_results,
          "results/05_Validation_Results.csv", row.names = FALSE)
write.csv(data.frame(Gene = available_genes),
          "results/05_Final_Gene_Signature.csv", row.names = FALSE)

cat("✓ Saved:\n")
cat("  - results/05_Validation_Final.rds\n")
cat("  - results/05_Validation_Results.csv\n")
cat("  - results/05_Final_Gene_Signature.csv\n\n")

# ============================================
# [5.10] Final summary
# ============================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║          ✓ STEP 5 COMPLETE — VALIDATION DONE                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("FINAL GENE SIGNATURE (%d genes):\n", length(available_genes)))
cat(paste(rep("─", 62), collapse = ""), "\n")
for (i in seq_along(available_genes)) {
  cat(sprintf("  %2d. %s\n", i, available_genes[i]))
}
cat("\n")

cat("VALIDATION ON INDEPENDENT DATASETS:\n")
cat(paste(rep("─", 62), collapse = ""), "\n")

for (i in seq_len(nrow(validation_results))) {
  cat(sprintf("  %s (N=%d, PE=%d, Control=%d):\n",
              validation_results$Dataset[i],
              validation_results$N_Samples[i],
              validation_results$N_PE[i],
              validation_results$N_Control[i]))
  cat(sprintf("    Genes found : %d / %d\n",
              validation_results$N_Genes[i], length(best_genes)))
  cat(sprintf("    AUC         : %s\n", validation_results$AUC[i]))
  cat(sprintf("    Sensitivity : %s\n", validation_results$Sensitivity[i]))
  cat(sprintf("    Specificity : %s\n\n", validation_results$Specificity[i]))
}

cat(paste(rep("═", 62), collapse = ""), "\n")
cat("✓ Ready for Figure generation!\n")
cat(paste(rep("═", 62), collapse = ""), "\n")

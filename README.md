# =============================================================================
# ENHANCED LATENT CLASS ANALYSIS (LCA) CLUSTERING for BINARY DATA
# =============================================================================
# Author : Enhanced Analysis Pipeline
# Created: 2025
# Purpose: Identify and characterize distinct subpopulations in CMML using 
#          Latent Class Analysis (LCA) based on baseline NGS (binary 0/1).
# 
# Key Enhancements:
# - LCA-based clustering (BayesLCA) for baseline mutations
# - Feature set: Baseline NGS (NGS_1)
# - Robust clinical and NGS feature engineering
# - Statistical outputs: BIC, class membership probabilities, defining features per class
# - Visualizations: BIC plots, class membership, density plots, heatmaps, proximity graph for cluster, charecters per cluster like: baseline blood counts, cytogenetics, subtype, and progression
# - Cluster summary table: clinical/molecular features by class
# - Survival analysis of LCA clusters
# - Silhouette plot for cluster validation (using binary distance)
# - Save all outputs to CSV and images
# - Comprehensive summary and recommendations at the end
# =============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# =============================================================================

library(BayesLCA)

required_packages <- c(
  "data.table", "BayesLCA", "ggplot2", "cluster", "pheatmap", "RColorBrewer", 
  "survival", "survminer", "dplyr", "reshape2", "gridExtra","readr", "factoextra"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

message("=== ENHANCED LCA CLUSTERING (BayesLCA) ===")

# =============================================================================
# SECTION 2: DATA LOADING AND PREPARATION
# =============================================================================

clinical_file <- "Sheet1.csv"
mutation_file <- "CMML_Project_2.csv"
df_clinical <- fread(clinical_file, stringsAsFactors = FALSE)
df_ngs <- fread(mutation_file, stringsAsFactors = FALSE)

# Ensure MRN is character for merging
df_clinical$MRN <- as.character(df_clinical$MRN)
df_ngs$MRN <- as.character(df_ngs$MRN)
makedata <- function(D){
  Df <- data.frame(D,check.names = FALSE)
  nr=nrow(Df);nc=ncol(Df)
  Df[Df==""]<- NA
  r <- apply(Df,1,function(x) sum(is.na(x)))
  c <- apply(Df,2,function(x) sum(is.na(x)))
  if(length(which(r==nc))){Df <- Df[-which(r==nc),]}
  if(length(which(c==nr))){Df <- Df[,-which(c==nr)]}
  return(Df)
}
df_clinical <- makedata(df_clinical)
df_ngs <- makedata(df_ngs)
d <- merge(df_clinical, df_ngs, by = "MRN", all.x = TRUE)
df <- data.table(d)

# =============================================================================
# SECTION 3: FEATURE ENGINEERING
# =============================================================================

clean_numeric <- function(x) {
  x <- gsub("[^0-9.\\-eE]", "", as.character(x))
  suppressWarnings(as.numeric(x))
}

# Use ONLY selected VAF columns (user provided)
selected_vaf_cols <- c(
  "NGS_1ABL1 VAF", "NGS_1ASXL1 VAF", "NGS_1BCOR VAF", "NGS_1BCORL1 VAF", "NGS_1BRAF VAF", "NGS_1CALR VAF", "NGS_1CBL VAF",
  "NGS_1CDKN2A VAF", "NGS_1CEBPA VAF", "NGS_1CSF3R VAF", "NGS_1CUX1 VAF", "NGS_1DDX41 VAF", "NGS_1DNMT3A VAF",
  "NGS_1EED VAF", "NGS_1ETNK1 VAF", "NGS_1ETV6 VAF", "NGS_1EZH2 VAF", "NGS_1FBXW7 VAF", "NGS_1FLT3 VAF", "NGS_1GNAS VAF",
  "NGS_1GATA2 VAF", "NGS_1G0S VAF", "NGS_1IDH1 VAF", "NGS_1IDH2 VAF", "NGS_1IKZF1 VAF", "NGS_1JAK2 VAF", "NGS_1JAK3 VAF",
  "NGS_1KDM6A VAF", "NGS_1KIT VAF", "NGS_1KMT2A VAF", "NGS_1KRAS VAF", "NGS_1LUC7L2 VAF", "NGS_1MPL VAF", "NGS_1MYD88 VAF",
  "NGS_1NF1 VAF", "NGS_1NOTCH1 VAF", "NGS_1NPM1 VAF", "NGS_1NRAS VAF", "NGS_1PHF6 VAF", "NGS_1PIGA VAF", "NGS_1PPM1D VAF",
  "NGS_1PRPF8 VAF", "NGS_1PTEN VAF", "NGS_1PTPN11 VAF", "NGS_1RAD21 VAF", "NGS_1RIT1 VAF", "NGS_1RUNX1 VAF",
  "NGS_1SETBP1 VAF", "NGS_1SF3B1 VAF", "NGS_1SH2B3 VAF", "NGS_1SMC1A VAF", "NGS_1SMC3 VAF", "NGS_1SRSF2 VAF",
  "NGS_1STAG2 VAF", "NGS_1STAT3 VAF", "NGS_1STAT5B VAF", "NGS_1SUZ12 VAF", "NGS_1TET2 VAF", "NGS_1TP53 VAF",
  "NGS_1U2AF1 VAF", "NGS_1U2AF2 VAF", "NGS_1WT1 VAF", "NGS_1ZRSR2 VAF"
)

# Ensure only columns present in your data are used
selected_vaf_cols <- intersect(selected_vaf_cols, names(df))
if (length(selected_vaf_cols) == 0) {
  stop("No selected VAF columns found in your input data. Check your column names!")
}

vaf_matrix <- df[, ..selected_vaf_cols]
vaf_matrix[] <- lapply(vaf_matrix, clean_numeric)

# Convert to binary: 0 if VAF==0, 1 if VAF>0 (NAs stay NA)
mutation_matrix <- as.data.frame(vaf_matrix)
mutation_matrix[] <- lapply(mutation_matrix, function(x) ifelse(is.na(x), NA, as.integer(x > 0)))

# --- KEY CHANGE: Keep only genes mutated in at least 5% of patients ---
n_patients <- nrow(mutation_matrix)
min_patients <- ceiling(0.04* n_patients)
keep_cols <- colSums(mutation_matrix == 1, na.rm = TRUE) >= min_patients
mutation_matrix <- mutation_matrix[, keep_cols, drop = FALSE]

# Keep only patients with at least 1 mutation
keep_rows <- rowSums(mutation_matrix == 1, na.rm = TRUE) >= 1
mutation_matrix <- mutation_matrix[keep_rows, , drop = FALSE]
df_for_lca <- df[keep_rows, ]

cat("Genes retained after filter (min_patients=5%):", sum(keep_cols), "\n")
cat("Patients retained after filter (min_patients=5%):", nrow(mutation_matrix), "\n")

# Remove columns with all NA or all 0
keep_cols <- colSums(mutation_matrix == 1, na.rm = TRUE) >= 1
mutation_matrix <- mutation_matrix[, keep_cols, drop = FALSE]

# Remove rows with all NA or all 0
keep_rows <- rowSums(mutation_matrix == 1, na.rm = TRUE) >= 1
mutation_matrix <- mutation_matrix[keep_rows, , drop = FALSE]

# Impute any remaining NAs as 0
mutation_matrix[is.na(mutation_matrix)] <- 0

# =============================================================================
# SECTION 4: LCA CLUSTERING WITH BayesLCA & STATISTICAL OUTPUTS
# =============================================================================

message("Fitting Latent Class Model (LCA) using BayesLCA for filtering (min_patients=5%)...")
max_classes <- 6
bic_vec <- rep(NA_real_, max_classes)
lca_fits <- vector("list", max_classes)
for (k in 3:max_classes) {
  set.seed(123)
  tryCatch({
    fit <- blca.em(as.matrix(mutation_matrix), G = k)
    bic_vec[k] <- fit$BIC
    lca_fits[[k]] <- fit
    cat(sprintf("Classes: %d, BIC: %.2f\n", k, fit$BIC))
  }, error=function(e) cat(sprintf("Classes: %d, ERROR: %s\n", k, e$message)))
}
bic_table <- data.frame(Classes = 3:max_classes, BIC = bic_vec[3:max_classes])
write.csv(bic_table, "LCA_BIC_table_minPatients5percent.csv", row.names = FALSE)
bic_plot <- ggplot(bic_table, aes(x = Classes, y = BIC)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Model Selection Criteria (BayesLCA, min_patients=5%)", x = "# of Classes", y = "BIC") +
  theme_minimal()
ggsave("LCA_BIC_plot_minPatients5percent.png", plot = bic_plot, width = 8, height = 5, dpi = 300)
best_k <- which.min(bic_vec)
best_fit <- lca_fits[[best_k]]

lca_clusters <- apply(best_fit$Z, 1, which.max)
df_for_lca$LCA_cluster <- factor(lca_clusters)
membership_probs <- as.data.frame(best_fit$Z)
colnames(membership_probs) <- paste0("Class", 1:best_k)
df_for_lca <- cbind(df_for_lca, membership_probs)
feature_names <- colnames(mutation_matrix)
class_means <- aggregate(mutation_matrix, by = list(Class = lca_clusters), mean, na.rm = TRUE)
write.csv(class_means, "LCA_class_feature_means_minPatients5percent.csv", row.names = FALSE)
message(paste("Optimal number of clusters selected by BIC (min_patients=5%):", best_k))
print(table(df_for_lca$LCA_cluster))

# Silhouette Plot
if(length(unique(lca_clusters)) > 1) {
  bin_dist <- stats::dist(mutation_matrix, method = "binary")  # FIXED: Explicit base R stats::dist
  sil <- silhouette(lca_clusters, bin_dist)
  avg_sil <- mean(sil[, 3])
  sil_plot <- factoextra::fviz_silhouette(sil)
  ggsave("LCA_silhouette_plot_minPatients5percent.png", plot = sil_plot, width = 8, height = 5, dpi = 300)
  if (!is.na(avg_sil)) {
    cat(sprintf("\nAverage silhouette width (min_patients=5%%): %.3f\n", avg_sil))
  } else {
    cat("\nAverage silhouette width (min_patients=2%%): NA\n")
  }
} else {
  message("Silhouette plot not generated: only one cluster detected.")
}

# =============================================================================
# SECTION 5: VISUALIZATION
# =============================================================================

# 5.1 Heatmap of feature means by class
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)
class_mean_mat <- as.matrix(class_means[, -1])
rownames(class_mean_mat) <- paste0("C", class_means$Class)
pheatmap(class_mean_mat,
         cluster_rows = FALSE, cluster_cols = TRUE,
         main = "Cluster Feature Means (LCA Clusters, min_patients=2%)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "LCA_feature_heatmap_minPatients5percent.png", width = 12, height = 6)

# 5.2 Heatmap of all patients: features x patient, annotated by cluster
patient_heatmap <- mutation_matrix
rownames(patient_heatmap) <- df_for_lca$MRN
annotation_row <- data.frame(LCA_cluster = df_for_lca$LCA_cluster)
rownames(annotation_row) <- rownames(patient_heatmap)
pheatmap(patient_heatmap,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_row = annotation_row,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "All Patients - Feature Heatmap by LCA Cluster (min_patients=2%)",
         filename = "LCA_allpatients_heatmap_minPatients5percent.png", width = 15, height = 8)

# 5.3 Centroid Distance Heatmap (Proximity of LCA Clusters)
centroid_matrix <- as.matrix(class_means[,-1])
rownames(centroid_matrix) <- paste0("C", class_means$Class)
centroid_dist <- dist(centroid_matrix)
centroid_dist_mat <- as.matrix(centroid_dist)
pheatmap(
  centroid_dist_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Centroid Distance Heatmap (LCA Clusters, min_patients=2%)",
  filename = "LCA_centroid_distance_heatmap_minPatients5percent.png",
  width = 8, height = 6
)

# =============================================================================
# SECTION 9: SUMMARY AND RECOMMENDATIONS
# =============================================================================
cat("\n=== SUMMARY ===\n")
cat("Analysis method: Latent Class Analysis (BayesLCA) on baseline mutational, clinical, and cytogenetic features. Filtering: genes mutated in >=2% of patients, patients with >=1 mutation.\n")
cat("Optimal clusters (BIC):", best_k, "\n")
cat("Patients per cluster:\n")
print(table(df_for_lca$LCA_cluster))
cat("\nModel fit statistics: see LCA_BIC_table_minPatients5percent.csv for details.\n")
cat("Defining feature means: see LCA_class_feature_means_minPatients5percent.csv.\n")
cat("\nClustered data saved to: CMML_clinical_LCA_clustered_minPatients5percent.csv.\n")
cat("\nSilhouette plot saved as LCA_silhouette_plot_minPatients5percent.png.\n")
cat("\nHeatmaps saved as LCA_feature_heatmap_minPatients5percent.png, LCA_allpatients_heatmap_minPatients5percent.png, and LCA_centroid_distance_heatmap_minPatients5percent.png.\n")

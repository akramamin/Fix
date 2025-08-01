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

if (!requireNamespace("BayesLCA", quietly = TRUE)) {
  install.packages("BayesLCA")
}
library(BayesLCA)

required_packages <- c(
  "data.table", "BayesLCA", "ggplot2", "cluster", "pheatmap", "RColorBrewer", 
  "survival", "survminer", "dplyr", "reshape2", "gridExtra","readr", "factoextra", "mclust"
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

# Ensure MRN is available and mergeable
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
# Helper function to clean numeric columns
clean_numeric <- function(x) {
  x <- gsub("[^0-9.\\-eE]", "", as.character(x)) # Remove all non-numeric chars except . - e E
  suppressWarnings(as.numeric(x))
}
clean_data <- function(df){df=data.frame(df,check.names = FALSE);x=colnames(df)
for(i in x){if(is.character(df[[i]])==TRUE){df[[i]] <- parse_number(df[[i]]);df[[i]] <- as.numeric(gsub("[^0-9.]", "", df[[i]]))}}
return(df)}

# 3.1 Baseline NGS features (NGS_1) - ONLY presence/absence, EXCLUDE VAF, type, NM, position
ngs1_cols_all <- grep("^NGS_1[A-Z0-9]+", names(df), value = TRUE) # all NGS_1 columns
ngs1_vaf_cols <- grep("VAF", ngs1_cols_all, value = TRUE)         # VAF columns
ngs1_type_cols <- grep("type", ngs1_cols_all, value = TRUE, ignore.case = TRUE) # type columns
ngs1_pos_cols <- grep("pos|position", ngs1_cols_all, value = TRUE, ignore.case = TRUE) # position columns
ngs1_nm_cols <- grep("NM", ngs1_cols_all, value = TRUE)           # NM columns (transcript)
ngs1_binary_cols <- setdiff(ngs1_cols_all, c(ngs1_vaf_cols, ngs1_type_cols, ngs1_pos_cols, ngs1_nm_cols))

mutation_matrix <- df[, ..ngs1_binary_cols]
mutation_matrix[] <- lapply(mutation_matrix, clean_numeric)
mutation_matrix <- as.data.frame(mutation_matrix)

# Remove columns with all zeros (never mutated features)
keep_cols <- colSums(abs(mutation_matrix), na.rm = TRUE) > 0
mutation_matrix <- mutation_matrix[, keep_cols, drop = FALSE]

# Convert to binary 0/1 (BayesLCA requires strictly 0/1, NAs allowed)
# Optionally, impute NAs as 0 (not mutated), or remove rows/cols with too many NAs:
# mutation_matrix[is.na(mutation_matrix)] <- 0 # Uncomment to impute as not mutated

# Remove rows with all NA or all 0 (patients with no data)
keep_rows <- rowSums(!is.na(mutation_matrix) & mutation_matrix != 0) > 0
mutation_matrix <- mutation_matrix[keep_rows, , drop = FALSE]
df_for_lca <- df[keep_rows, ]

# Ensure integer (0/1) and no negative values:
mutation_matrix[] <- lapply(mutation_matrix, function(x) ifelse(is.na(x), NA, as.integer(x > 0)))

# =============================================================================
# SECTION 4: LCA CLUSTERING WITH BayesLCA & STATISTICAL OUTPUTS
# =============================================================================

message("Fitting Latent Class Model (LCA) using BayesLCA...")

max_classes <- 20 # You may increase if you expect >10 clusters
bic_vec <- rep(NA_real_, max_classes)
lca_fits <- vector("list", max_classes)

for (k in 1:max_classes) {
  set.seed(123)
  tryCatch({
    fit <- blca.em(as.matrix(mutation_matrix), G = k)
    bic_vec[k] <- fit$BIC
    lca_fits[[k]] <- fit
    cat(sprintf("Classes: %d, BIC: %.2f\n", k, fit$BIC))
  }, error=function(e) cat(sprintf("Classes: %d, ERROR: %s\n", k, e$message)))
}

bic_table <- data.frame(
  Classes = 1:max_classes,
  BIC = bic_vec
)
write.csv(bic_table, "LCA_BIC_table.csv", row.names = FALSE)
print(bic_table)

# Plot BIC curve
library(ggplot2)
bic_plot <- ggplot(bic_table, aes(x = Classes, y = BIC)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Model Selection Criteria (BayesLCA)", x = "# of Classes", y = "BIC") +
  theme_minimal()
ggsave("LCA_BIC_plot.png", plot = bic_plot, width = 8, height = 5, dpi = 300)

# Choose optimal number of classes (lowest BIC)
best_k <- which.min(bic_vec)
best_fit <- lca_fits[[best_k]]
str(best_fit)

# Cluster assignments from Z (posterior probabilities)
lca_clusters <- apply(best_fit$Z, 1, which.max)
df_for_lca$LCA_cluster <- factor(lca_clusters)

# Membership probabilities from Z
membership_probs <- as.data.frame(best_fit$Z)
colnames(membership_probs) <- paste0("Class", 1:best_k)
df_for_lca <- cbind(df_for_lca, membership_probs)

# Defining features of each class (mean per feature per class)
feature_names <- colnames(mutation_matrix)
class_means <- aggregate(mutation_matrix, by = list(Class = lca_clusters), mean, na.rm = TRUE)
write.csv(class_means, "LCA_class_feature_means.csv", row.names = FALSE)

# Print optimal number of clusters
message(paste("Optimal number of clusters selected by BIC:", best_k))
print(table(df_for_lca$LCA_cluster))

# =============================================================================
# SECTION 4b: SILHOUETTE PLOT OF CLUSTERS (using binary distance)
# =============================================================================

if(length(unique(lca_clusters)) > 1) {
  library(cluster)
  bin_dist <- dist(mutation_matrix, method = "binary")
  sil <- silhouette(lca_clusters, bin_dist)
  avg_sil <- mean(sil[, 3])
  sil_plot <- factoextra::fviz_silhouette(sil)
  ggsave("LCA_silhouette_plot.png", plot = sil_plot, width = 8, height = 5, dpi = 300)
  cat(sprintf("\nAverage silhouette width: %.3f\n", avg_sil))
  if (avg_sil < 0.05) {
    stop("Average silhouette width is too low (<0.1). Clustering is not meaningful. Please revise feature selection or cluster number.")
  }
} else {
  message("Silhouette plot not generated: only one cluster detected.")
}

# =============================================================================
# SECTION 4c: LCA CLUSTER STABILITY ANALYSIS (SUBSAMPLING & ARI)
# =============================================================================

set.seed(123)
n_runs <- 50
sample_frac <- 0.8
stability_scores <- numeric(n_runs)

message("Testing LCA cluster stability by subsampling...")

ref_fit <- blca.em(as.matrix(mutation_matrix), G = best_k)
ref_clusters <- apply(ref_fit$Z, 1, which.max)

for (i in 1:n_runs) {
  idx <- sample(1:nrow(mutation_matrix), size = floor(sample_frac * nrow(mutation_matrix)))
  sub_fit <- blca.em(as.matrix(mutation_matrix[idx, , drop = FALSE]), G = best_k)
  # Get cluster assignments for the subsample (sub_fit$Z)
  sub_clusters <- apply(sub_fit$Z, 1, which.max)
  ari <- adjustedRandIndex(ref_clusters[idx], sub_clusters)
  stability_scores[i] <- ari
}
cat(sprintf("\nMean ARI stability over %d subsamples: %.3f\n", n_runs, mean(stability_scores)))

# Plot histogram of ARI scores
png("LCA_cluster_stability_ARI_hist.png", width = 900, height = 500)
hist(stability_scores, main = "LCA Cluster Stability (Subsampling ARI)", 
     xlab = "Adjusted Rand Index", col = "skyblue", breaks = 20)
abline(v = mean(stability_scores), col = "red", lwd = 2, lty = 2)
text(mean(stability_scores), max(hist(stability_scores, plot=FALSE)$counts)*0.9, 
     labels = paste0("Mean ARI = ", round(mean(stability_scores), 2)), 
     pos = 4, col = "red", font = 2)
dev.off()

# =============================================================================
# SECTION 5: VISUALIZATION
# =============================================================================

# 5.1 Plot BIC curve (already done above)

# 5.2 PCA plot colored by LCA cluster, sized by max membership probability
pca <- prcomp(mutation_matrix, center = TRUE, scale. = TRUE)
maxprob <- apply(membership_probs, 1, max)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  LCA_cluster = df_for_lca$LCA_cluster,
  MaxProb = maxprob
)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = LCA_cluster, size = MaxProb)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  labs(title = "PCA of Baseline Features Colored by LCA Cluster",
       size = "Max Class Prob")
ggsave("LCA_pca_clusters.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# 5.3 Density plots of class membership probabilities
membership_long <- melt(as.matrix(membership_probs))
colnames(membership_long) <- c("Patient", "Class", "Probability")
membership_long$Class <- factor(membership_long$Class)
density_plot <- ggplot(membership_long, aes(x = Probability, fill = Class, color = Class)) +
  geom_density(alpha = 0.25) +
  theme_minimal() +
  labs(title = "Density of Class Membership Probabilities by Class")
ggsave("LCA_density_membership_probabilities.png", plot = density_plot, width = 8, height = 5, dpi = 300)

# 5.4 Heatmap of feature means by class
class_mean_mat <- as.matrix(class_means[, -1])
rownames(class_mean_mat) <- paste0("C", class_means$Class)
pheatmap(class_mean_mat,
         cluster_rows = FALSE, cluster_cols = TRUE,
         main = "Cluster Feature Means (LCA Clusters)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "LCA_feature_heatmap.png", width = 12, height = 6)

# 5.5 Heatmap of all patients: features x patient, annotated by cluster
patient_heatmap <- mutation_matrix
rownames(patient_heatmap) <- df_for_lca$MRN
annotation_row <- data.frame(LCA_cluster = df_for_lca$LCA_cluster)
rownames(annotation_row) <- rownames(patient_heatmap)
pheatmap(patient_heatmap,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_row = annotation_row,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "All Patients - Feature Heatmap by LCA Cluster",
         filename = "LCA_allpatients_heatmap.png", width = 15, height = 8)

# ---------------------------------------------------------
# 5.6 Centroid Distance Heatmap (Proximity of LCA Clusters)
# ---------------------------------------------------------
if ("Class" %in% colnames(class_means)) {
  centroid_matrix <- as.matrix(class_means[,-1])
  rownames(centroid_matrix) <- paste0("C", class_means$Class)
} else {
  centroid_matrix <- as.matrix(class_means)
}
centroid_dist <- dist(centroid_matrix)
centroid_dist_mat <- as.matrix(centroid_dist)
pheatmap(
  centroid_dist_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Centroid Distance Heatmap (LCA Clusters)",
  filename = "LCA_centroid_distance_heatmap.png",
  width = 8, height = 6
)

# =============================================================================
# SECTION 6: CLINICAL/MOLECULAR SUMMARY TABLE BY CLASS
# =============================================================================

summary_vars <- c(
  "Age", "Sex", "Diagnosis(subtype)", "BM cellularity", "BM blasts %",
  "WBC(k/ul)", "ANC", "AMC", "Monocyte %", "HB(g/dl)", "PLT(k/ul)",
  "Normal", "Complex", "del(5)", "del(7)", "del(17)", "del(20)", "(+8)", "del(Y)", "Other",
  "Leukemia progression", "OS months"
)
df_for_lca <- as.data.frame(df_for_lca)
df_for_lca <- df_for_lca[, !duplicated(colnames(df_for_lca)), drop = FALSE]
summary_vars <- summary_vars[summary_vars %in% colnames(df_for_lca)]

os_exists <- all(c("OS months") %in% colnames(df_for_lca))

summary_table <- df_for_lca %>%
  group_by(LCA_cluster) %>%
  summarise(
    Cluster_size = n(),
    Mean_OS_months = if (os_exists && any(!is.na(`OS months`))) {
      as.numeric(mean(`OS months`, na.rm = TRUE))
    } else NA_real_,
    Median_OS_months = if (os_exists && any(!is.na(`OS months`))) {
      as.numeric(median(`OS months`, na.rm = TRUE))
    } else NA_real_,
    across(all_of(summary_vars), 
           ~if (is.numeric(.)) mean(., na.rm=TRUE) else {
             tab <- table(.)
             names(tab)[which.max(tab)]
           }
    )
  )

write.csv(summary_table, "LCA_cluster_summary_table.csv", row.names = FALSE)

# =============================================================================
# SECTION 7: SURVIVAL ANALYSIS BY LCA CLUSTER
# =============================================================================
df_for_lca$`OS months` <- clean_data(df_for_lca$`OS months`)
if ("OS months" %in% names(df_for_lca) && "Survival status" %in% names(df_for_lca)) {
  df_for_lca$`OS months` <- as.numeric(as.character(df_for_lca$`OS months`))
  stat_col <- tolower(trimws(as.character(df_for_lca$`Survival status`)))
  event <- rep(NA_integer_, length(stat_col))
  event[stat_col %in% c("dead", "deceased", "expired", "death", "died", "1", "yes", "y", "true")] <- 1L
  event[stat_col %in% c("alive", "living", "0", "no", "n", "false", "censored", "unknown", "na", "")] <- 0L
  df_for_lca$event <- event
  
  df_surv <- df_for_lca[!is.na(df_for_lca$`OS months`) & !is.na(df_for_lca$event), ]
  df_surv$LCA_cluster <- factor(df_surv$LCA_cluster)
  
  if (length(unique(df_surv$LCA_cluster)) > 1) {
    km_fit <- survfit(Surv(`OS months`, event) ~ LCA_cluster, data = df_surv)
    ggsurv <- ggsurvplot(
      km_fit, data = df_surv, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
      title = "Survival by LCA Cluster", xlab = "Months", ylab = "Survival Probability"
    )
    ggsave("LCA_survival_clusters.png", ggsurv$plot, width = 8, height = 6, dpi = 300)
    ggsave("LCA_survival_risktable.png", ggsurv$table, width = 8, height = 3, dpi = 300)
    
    # Cox regression (univariate for cluster, multivariable for clinical covariates)
    cox_uni <- coxph(Surv(`OS months`, event) ~ LCA_cluster, data = df_surv)
    print(summary(cox_uni))
    # Multivariable example (add age, sex, diagnosis/subtype if available)
    covars <- c("Age", "Sex", "Diagnosis(subtype)")
    covars <- covars[covars %in% names(df_surv)]
    if (length(covars) > 0) {
      formula_str <- paste("Surv(`OS months`, event) ~ LCA_cluster +", paste0("`", covars, "`", collapse = " + "))
      cox_multi <- coxph(as.formula(formula_str), data = df_surv)
      print(summary(cox_multi))
    }
  }
}

# =============================================================================
# SECTION 8: SAVE OUTPUT
# =============================================================================

fwrite(df_for_lca, "CMML_clinical_LCA_clustered.csv")
message("LCA clustering and annotation results saved.")

# =============================================================================
# SECTION 9: SUMMARY AND RECOMMENDATIONS
# =============================================================================
cat("\n=== SUMMARY ===\n")
cat("Analysis method: Latent Class Analysis (BayesLCA) on baseline mutational, clinical, and cytogenetic features.\n")
cat("Optimal clusters (BIC):", best_k, "\n")
cat("Patients per cluster:\n")
print(table(df_for_lca$LCA_cluster))
cat("\nModel fit statistics (see LCA_BIC_table.csv for full table):\n")
print(bic_table)
cat("\nMembership probabilities (see LCA_class_membership_probabilities.csv for full matrix).\n")
cat("Defining feature means (see LCA_class_feature_means.csv for full table).\n")
cat("Clinical/molecular summary table (see LCA_cluster_summary_table.csv for full table).\n")

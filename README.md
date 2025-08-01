# =============================================================================
# ENHANCED FINITE MIXTURE MODEL (FMM) CLUSTERING
# =============================================================================
# Author : Enhanced Analysis Pipeline
# Created: 2025
# Purpose: Identify and characterize distinct subpopulations in CMML using 
#          Finite Mixture Models (FMM) based on baseline NGS, clinical, and 
#          cytogenetics, and key clinical outcomes. (NO trajectory analysis)
# 
# Key Enhancements:
# - FMM-based clustering (mclust) for baseline mutations
# - Feature set: Baseline NGS (NGS_1)
# - Robust clinical and NGS feature engineering
# - Statistical outputs: BIC/AIC/likelihood, class membership probabilities, defining features per class
# - Visualizations: BIC/AIC plots, class membership, density plots, heatmaps, proximity graph for cluster, charecters per cluster like: baseline blood counts, cytogenetics, subtype, and progression
# - Cluster summary table: clinical/molecular features by class
# - Survival analysis of FMM clusters
# - Silhouette plot for cluster validation
# - Save all outputs to CSV and images
# - Comprehensive summary and recommendations at the end
# =============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# =============================================================================

required_packages <- c(
  "data.table", "mclust", "ggplot2", "cluster", "pheatmap", "RColorBrewer", 
  "survival", "survminer", "dplyr", "reshape2", "gridExtra","readr", "factoextra"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

message("=== ENHANCED FMM CLUSTERING ===")

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
# DO NOT IMPUTE NA for mutation_matrix!
# 1 = mutated, 0 = wildtype, NA = not tested
mutation_matrix <- as.data.frame(mutation_matrix)

# 3.2 FMM input (add clinical/cyto features here if desired)
fmm_input <- mutation_matrix

# Remove columns with all zeros BEFORE scaling
keep_cols <- colSums(abs(fmm_input), na.rm = TRUE) > 0
fmm_input <- fmm_input[, keep_cols, drop = FALSE]

# Scale all features (mean 0, sd 1)
fmm_input <- scale(fmm_input)
library(caret)
nzv <- nearZeroVar(fmm_input)
if (length(nzv) > 0) {
  fmm_input <- fmm_input[, -nzv, drop = FALSE]
}
# Remove rows with all zeros (patients with no data)
keep_rows <- rowSums(is.na(fmm_input) | fmm_input == 0, na.rm = TRUE) < ncol(fmm_input)
fmm_input <- fmm_input[keep_rows, ]
df_for_fmm <- df[keep_rows, ]

# After scaling, check for NA/NaN/Inf
if (any(!is.finite(as.matrix(fmm_input)))) {
  message("Non-finite values present after scaling, replacing with 0")
  fmm_input[!is.finite(fmm_input)] <- 0
}

# =============================================================================
# SECTION 4: FMM CLUSTERING WITH MCLUST & STATISTICAL OUTPUTS
# =============================================================================

message("Fitting Finite Mixture Model (FMM) using mclust...")

fmm_fit <- Mclust(fmm_input, G = 1:13)

fmm_clusters <- fmm_fit$classification
df_for_fmm$FMM_cluster <- factor(fmm_clusters)

# 4.1 Statistical outputs: BIC, AIC, log-likelihood, n_parameters
bic_values <- fmm_fit$BIC
loglik_values <- fmm_fit$loglik
n_parameters <- nMclustParams(fmm_fit$modelName, fmm_fit$d, fmm_fit$G)
aic_values <- 2 * n_parameters - 2 * loglik_values

# Get matrix of values
bic_matrix <- as.matrix(bic_values)  # Just the numeric part

# Get the row names (number of clusters G) and column names (models)
G_values <- as.numeric(rownames(bic_matrix))
model_names <- colnames(bic_matrix)

# Convert to long format manually
bic_df <- data.frame(
  G = rep(G_values, each = length(model_names)),
  Model = rep(model_names, each = length(G_values)),
  BIC = as.vector(bic_matrix)
)
bic_dff=data.frame(1:length(G_values))
for (i in unique(bic_df$Model)) {
  bic_dff=data.frame(bic_dff,bic_df$BIC[which(bic_df$Model==i)])
}
colnames(bic_dff) <- c("G",unique(bic_df$Model))

# Save BIC/AIC/loglik values per G
bic_aic_loglik <- data.frame(
  G = 1:13,
  BIC=bic_values,
  AIC = aic_values,
  LogLikelihood = loglik_values 
)
bic_aic_logli=list()
bic_aic_logli$BIC=bic_dff
bic_aic_logli$AIC = aic_values
bic_aic_logli$LogLikelihood = loglik_values
bic_aic_logli$G=1:13

write.csv(bic_aic_logli, "FMM_BIC_AIC_LogLik.csv", row.names = FALSE)

# 4.2 Class membership probabilities for each patient
membership_probs <- fmm_fit$z
df_for_fmm <- cbind(df_for_fmm, membership_probs)
write.csv(df_for_fmm[, c("MRN", "FMM_cluster", colnames(membership_probs))], "FMM_class_membership_probabilities.csv", row.names = FALSE)

# 4.3 Defining features of each class (mean per feature per class)
feature_names <- colnames(fmm_input)
class_means <- aggregate(fmm_input, by = list(Class = fmm_clusters), mean)
write.csv(class_means, "FMM_class_feature_means.csv", row.names = FALSE)

# Print optimal number of clusters and class sizes
message(paste("Optimal number of clusters selected by BIC:", fmm_fit$G))
print(table(df_for_fmm$FMM_cluster))

# =============================================================================
# SECTION 4b: SILHOUETTE PLOT OF CLUSTERS
# =============================================================================

if(length(unique(fmm_clusters)) > 1) {
  sil <- silhouette(fmm_clusters, dist(fmm_input))
  avg_sil <- mean(sil[, 3])
  sil_plot <- factoextra::fviz_silhouette(sil)
  ggsave("FMM_silhouette_plot.png", plot = sil_plot, width = 8, height = 5, dpi = 300)
  cat(sprintf("\nAverage silhouette width: %.3f\n", avg_sil))
  if (avg_sil < 0.1) {
    stop("Average silhouette width is too low (<0.1). Clustering is not meaningful. Please revise feature selection or cluster number.")
  }
} else {
  message("Silhouette plot not generated: only one cluster detected.")
}

# =============================================================================
# SECTION 4c: FMM CLUSTER STABILITY ANALYSIS (SUBSAMPLING & ARI)
# =============================================================================

if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
library(mclust)

set.seed(123)
n_runs <- 100
sample_frac <- 0.8
stability_scores <- numeric(n_runs)

message("Testing FMM cluster stability by subsampling...")

# Use the optimal number of clusters (fmm_fit$G)
ref_fit <- Mclust(fmm_input, G = fmm_fit$G)
ref_clusters <- ref_fit$classification

for (i in 1:n_runs) {
  idx <- sample(1:nrow(fmm_input), size = floor(sample_frac * nrow(fmm_input)))
  sub_fit <- Mclust(fmm_input[idx, , drop = FALSE], G = fmm_fit$G)
  # ARI between original clusters and subsample for those patients
  ari <- adjustedRandIndex(ref_clusters[idx], sub_fit$classification)
  stability_scores[i] <- ari
}
cat(sprintf("\nMean ARI stability over %d subsamples: %.3f\n", n_runs, mean(stability_scores)))

# Plot histogram of ARI scores
png("FMM_cluster_stability_ARI_hist.png", width = 900, height = 500)
hist(stability_scores, main = "FMM Cluster Stability (Subsampling ARI)", 
     xlab = "Adjusted Rand Index", col = "skyblue", breaks = 20)
abline(v = mean(stability_scores), col = "red", lwd = 2, lty = 2)
text(mean(stability_scores), max(hist(stability_scores, plot=FALSE)$counts)*0.9, 
     labels = paste0("Mean ARI = ", round(mean(stability_scores), 2)), 
     pos = 4, col = "red", font = 2)
dev.off()

# =============================================================================
# SECTION 5: VISUALIZATION
# =============================================================================

# 5.1 Plot BIC and AIC for different G
bic_aic_df <- data.frame(bic_aic_logli)
library(ggplot2)
bic_plot <- ggplot(bic_aic_df, aes(x = G)) +
  geom_line(aes(y = BIC.EII, color = "BIC")) +
  geom_line(aes(y = AIC, color = "AIC")) +
  geom_point(aes(y = BIC.EII, color = "BIC")) +
  geom_point(aes(y = AIC, color = "AIC")) +
  scale_color_manual(values = c("BIC" = "blue", "AIC" = "darkred")) +
  labs(title = "Model Selection Criteria", y = "Value", x = "# of Classes") +
  theme_minimal()
ggsave("FMM_BIC_AIC_plot.png", plot = bic_plot, width = 8, height = 5, dpi = 300)

# 5.2 PCA plot colored by FMM cluster, sized by max membership probability
pca <- prcomp(fmm_input, center = TRUE, scale. = TRUE)
maxprob <- apply(membership_probs, 1, max)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  FMM_cluster = df_for_fmm$FMM_cluster,
  MaxProb = maxprob
)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = FMM_cluster, size = MaxProb)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  labs(title = "PCA of Baseline Features Colored by FMM Cluster",
       size = "Max Class Prob")
ggsave("FMM_pca_clusters.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# 5.3 Density plots of class membership probabilities
membership_long <- melt(membership_probs)
colnames(membership_long) <- c("Patient", "Class", "Probability")
membership_long$Class <- factor(membership_long$Class)
density_plot <- ggplot(membership_long, aes(x = Probability, fill = Class, color = Class)) +
  geom_density(alpha = 0.25) +
  theme_minimal() +
  labs(title = "Density of Class Membership Probabilities by Class")
ggsave("FMM_density_membership_probabilities.png", plot = density_plot, width = 8, height = 5, dpi = 300)

# 5.4 Heatmap of feature means by class
class_mean_mat <- as.matrix(class_means[, -1])
rownames(class_mean_mat) <- paste0("C", class_means$Class)
pheatmap(class_mean_mat,
         cluster_rows = FALSE, cluster_cols = TRUE,
         main = "Cluster Feature Means (FMM Clusters)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "FMM_feature_heatmap.png", width = 12, height = 6)

# 5.5 Heatmap of all patients: features x patient, annotated by cluster
patient_heatmap <- fmm_input
rownames(patient_heatmap) <- df_for_fmm$MRN
annotation_row <- data.frame(FMM_cluster = df_for_fmm$FMM_cluster)
rownames(annotation_row) <- rownames(patient_heatmap)
pheatmap(patient_heatmap,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_row = annotation_row,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "All Patients - Feature Heatmap by FMM Cluster",
         filename = "FMM_allpatients_heatmap.png", width = 15, height = 8)

# =============================================================================
# SECTION 6: CLINICAL/MOLECULAR SUMMARY TABLE BY CLASS
# =============================================================================

summary_vars <- c(
  "Age", "Sex", "Diagnosis(subtype)", "BM cellularity", "BM blasts %",
  "WBC(k/ul)", "ANC", "AMC", "Monocyte %", "HB(g/dl)", "PLT(k/ul)",
  "Normal", "Complex", "del(5)", "del(7)", "del(17)", "del(20)", "(+8)", "del(Y)", "Other",
  "Leukemia progression"
)
df_for_fmm <- as.data.frame(df_for_fmm) # Force to data.frame
df_for_fmm <- df_for_fmm[, !duplicated(colnames(df_for_fmm)), drop = FALSE]
summary_vars <- summary_vars[summary_vars %in% colnames(df_for_fmm)]
summary_table <- df_for_fmm %>%
  group_by(FMM_cluster) %>%
  summarise(across(all_of(summary_vars), ~if (is.numeric(.)) mean(., na.rm=TRUE) else {
    tab <- table(.)
    names(tab)[which.max(tab)]
  }))
write.csv(summary_table, "FMM_cluster_summary_table.csv", row.names = FALSE)

# =============================================================================
# SECTION 7: SURVIVAL ANALYSIS BY FMM CLUSTER
# =============================================================================
df_for_fmm$`OS months` <- clean_data(df_for_fmm$`OS months`)
if ("OS months" %in% names(df_for_fmm) && "Survival status" %in% names(df_for_fmm)) {
  df_for_fmm$`OS months` <- as.numeric(as.character(df_for_fmm$`OS months`))
  stat_col <- tolower(trimws(as.character(df_for_fmm$`Survival status`)))
  event <- rep(NA_integer_, length(stat_col))
  event[stat_col %in% c("dead", "deceased", "expired", "death", "died", "1", "yes", "y", "true")] <- 1L
  event[stat_col %in% c("alive", "living", "0", "no", "n", "false", "censored", "unknown", "na", "")] <- 0L
  df_for_fmm$event <- event
  
  df_surv <- df_for_fmm[!is.na(df_for_fmm$`OS months`) & !is.na(df_for_fmm$event), ]
  df_surv$FMM_cluster <- factor(df_surv$FMM_cluster)
  
  if (length(unique(df_surv$FMM_cluster)) > 1) {
    km_fit <- survfit(Surv(`OS months`, event) ~ FMM_cluster, data = df_surv)
    ggsurv <- ggsurvplot(
      km_fit, data = df_surv, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
      title = "Survival by FMM Cluster", xlab = "Months", ylab = "Survival Probability"
    )
    ggsave("FMM_survival_clusters.png", ggsurv$plot, width = 8, height = 6, dpi = 300)
    ggsave("FMM_survival_risktable.png", ggsurv$table, width = 8, height = 3, dpi = 300)
    
    # Cox regression (univariate for cluster, multivariable for clinical covariates)
    cox_uni <- coxph(Surv(`OS months`, event) ~ FMM_cluster, data = df_surv)
    print(summary(cox_uni))
    # Multivariable example (add age, sex, diagnosis/subtype if available)
    covars <- c("Age", "Sex", "Diagnosis(subtype)")
    covars <- covars[covars %in% names(df_surv)]
    if (length(covars) > 0) {
      formula_str <- paste("Surv(`OS months`, event) ~ FMM_cluster +", paste0("`", covars, "`", collapse = " + "))
      cox_multi <- coxph(as.formula(formula_str), data = df_surv)
      print(summary(cox_multi))
    }
  }
}

# =============================================================================
# SECTION 8: SAVE OUTPUT
# =============================================================================

fwrite(df_for_fmm, "CMML_clinical_FMM_clustered.csv")
message("FMM clustering and annotation results saved.")

# =============================================================================
# SECTION 9: SUMMARY AND RECOMMENDATIONS
# =============================================================================
cat("\n=== SUMMARY ===\n")
cat("Analysis method: FMM clustering (mclust) on baseline mutational, clinical, and cytogenetic features.\n")
cat("Optimal clusters (BIC):", fmm_fit$G, "\n")
cat("Patients per cluster:\n")
print(table(df_for_fmm$FMM_cluster))
cat("\nModel fit statistics (see FMM_BIC_AIC_LogLik.csv for full table):\n")
print(bic_aic_logli)
cat("\nMembership probabilities (see FMM_class_membership_probabilities.csv for full matrix).\n")
cat("Defining feature means (see FMM_class_feature_means.csv for full table).\n")
cat("Clinical/molecular summary table (see FMM_cluster_summary_table.csv for full table).\n")

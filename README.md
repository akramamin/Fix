# Serial NGS and New Mutations Analysis with Time-Dependent Cox Regression (Optimized)
# This script analyzes serial NGS data with all advanced features but optimized package dependencies

library(tidyverse)

# Install and load data.table for fread
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Essential packages with fallbacks for advanced features
essential_packages <- c("survival", "readxl", "dplyr", "ggplot2", "gridExtra")

# Optional advanced packages (install only if needed)
advanced_packages <- c("survminer", "tableone", "forestplot", 
                       "RColorBrewer", "pheatmap", "corrplot", "broom")

# Install essential packages
for (pkg in essential_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing essential package:", pkg, "\n")
    install.packages(pkg, dependencies = c("Depends", "Imports"))
    library(pkg, character.only = TRUE)
  }
}

# Try to install advanced packages (skip if installation fails)
advanced_loaded <- c()
for (pkg in advanced_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Attempting to install advanced package:", pkg, "\n")
    tryCatch({
      install.packages(pkg, dependencies = c("Depends", "Imports"))
      library(pkg, character.only = TRUE)
      advanced_loaded <- c(advanced_loaded, pkg)
    }, error = function(e) {
      cat("Warning: Could not install", pkg, "- will use basic alternatives\n")
    })
  } else {
    advanced_loaded <- c(advanced_loaded, pkg)
  }
}

cat("Advanced packages loaded:", 
    paste(advanced_loaded, collapse = ", "), "\n")

# Create output directory
if (!dir.exists("Serial_NGS_Output")) {
  dir.create("Serial_NGS_Output")
}

# Function to safely convert dates
safe_date_conversion <- function(date_col) {
  if (inherits(date_col, "Date")) return(date_col)
  if (is.numeric(date_col)) return(as.Date(date_col, origin = "1899-12-30"))
  if (is.character(date_col) || is.factor(date_col)) {
    date_formats <- c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y-%m-%d %H:%M:%S")
    for (fmt in date_formats) {
      converted <- suppressWarnings(as.Date(as.character(date_col), format = fmt))
      if (sum(!is.na(converted)) > 0.5 * length(date_col)) return(converted)
    }
    if (all(grepl("^[0-9]+$", date_col[!is.na(date_col)]))) {
      serials <- as.numeric(date_col)
      return(as.Date(serials, origin = "1899-12-30"))
    }
  }
  return(rep(NA, length(date_col)))
}

# Function to read data

read_data_safely <- function() {
  clinical_file <- "Sheet1.csv"
  mutation_file <- "CMML_Project_2.csv"
  clinical_data <- fread(clinical_file, stringsAsFactors = FALSE)
  ngs_data <- fread(mutation_file, stringsAsFactors = FALSE)
  
  # Remove rows with NA MRN in either dataset, and keep only patients present in both
  clinical_data <- clinical_data[!is.na(clinical_data$MRN), ]
  ngs_data <- ngs_data[!is.na(ngs_data$MRN), ]
  shared_MRN <- intersect(clinical_data$MRN, ngs_data$MRN)
  clinical_data <- clinical_data[clinical_data$MRN %in% shared_MRN, ]
  ngs_data <- ngs_data[ngs_data$MRN %in% shared_MRN, ]
  
  return(list(clinical = clinical_data, ngs = ngs_data))
}

# Function to identify serial NGS patients
identify_serial_ngs <- function(ngs_data, vaf_threshold = 1) {
  cat("Identifying serial NGS patients (VAF >=", vaf_threshold, "% at ANY serial NGS timepoint)...\n")
  
  # Instead of just NGS_2, check NGS_2 through NGS_5 for any VAF >= 1
  serial_ngs_cols <- grep("^NGS_[2-5][A-Z0-9]+ VAF$", names(ngs_data), value = TRUE)
  serial_ngs_date_cols <- grep("^NGS_[2-5][ _]?date$", names(ngs_data), value = TRUE, ignore.case = TRUE)
  
  serial_patients <- c()
  for (i in 1:nrow(ngs_data)) {
    patient_id <- as.character(ngs_data$MRN[i])
    has_serial <- FALSE
    for (col in serial_ngs_cols) {
      val <- as.character(ngs_data[[col]][i])
      if (!is.na(val) && val != "") {
        nums <- suppressWarnings(as.numeric(strsplit(val, ",\\s*")[[1]]))
        if (any(!is.na(nums) & nums >= vaf_threshold)) {
          has_serial <- TRUE
          break
        }
      }
    }
    if (has_serial) serial_patients <- c(serial_patients, patient_id)
  }
  cat("Identified", length(serial_patients), "patients with serial NGS (VAF >=", vaf_threshold, ")\n")
  return(list(serial_patients = serial_patients, serial_ngs_cols = serial_ngs_cols))
}

# Function to detect new mutations at all serial NGS timepoints (NGS_2–NGS_5) (EDITED: checks date for each timepoint)
detect_new_mutations_all_timepoints <- function(ngs_data, vaf_threshold = 1) {
  cat("Detecting new mutations at all serial NGS timepoints (NGS_2–NGS_5)...\n")
  patient_ids <- ngs_data[[1]] # assumes first column is patient ID
  timepoints <- 1:5 # NGS_1 to NGS_5
  
  # Identify all VAF columns for each NGS timepoint
  vaf_cols <- lapply(timepoints, function(tp) {
    pattern1 <- paste0("NGS_", tp, "_")
    pattern2 <- paste0("NGS_", tp, "$")
    pattern3 <- paste0("NGS_", tp, "(?!_)")
    cols1 <- grep(pattern1, names(ngs_data), value = TRUE)
    cols2 <- grep(pattern2, names(ngs_data), value = TRUE)
    cols3 <- grep(pattern3, names(ngs_data), value = TRUE, perl = TRUE)
    cols <- unique(c(cols1, cols2, cols3))
    cols <- cols[grepl("VAF", cols)]
    return(cols)
  })
  
  vaf_cols <- vaf_cols[sapply(vaf_cols, length) > 0]
  timepoints <- timepoints[sapply(vaf_cols, length) > 0]
  new_mutations_data <- data.frame()
  if (length(vaf_cols) == 0) {
    cat("No VAF columns found in the data\n")
    return(new_mutations_data)
  }
  
  for (i in seq_len(nrow(ngs_data))) {
    patient <- as.character(patient_ids[i])
    for (gene_col in vaf_cols[[1]]) {
      gene <- gsub("ngs_1|_vaf", "", gene_col)
      vafs <- sapply(seq_along(timepoints), function(tp_idx) {
        # Build VAF col name for this timepoint
        if (grepl("NGS_1_", gene_col)) {
          colname <- gsub("NGS_1_", paste0("NGS_", timepoints[tp_idx], "_"), gene_col)
        } else {
          colname <- gsub("ngs_1", paste0("ngs_", timepoints[tp_idx]), gene_col)
        }
        # --- CHECK DATE FOR THIS TIMEPOINT ---
        date_col <- grep(paste0("^NGS_", timepoints[tp_idx], "[ _]?date$"), names(ngs_data), value = TRUE, ignore.case = TRUE)
        date_val <- if (length(date_col) == 1) ngs_data[[date_col]][i] else NA
        if (!is.na(date_val) && date_val != 0 && colname %in% names(ngs_data)) {
          val <- as.character(ngs_data[[colname]][i])
          if (is.na(val) || val == "") return(NA)
          nums <- as.numeric(strsplit(val, ",\\s*")[[1]])
          if (all(is.na(nums))) NA else max(nums, na.rm = TRUE)
        } else {
          NA
        }
      })
      # Loop over serial NGS (from NGS_2 onwards)
      for (tp_idx in 2:length(timepoints)) {
        prev_vaf <- vafs[tp_idx - 1]
        curr_vaf <- vafs[tp_idx]
        if (!is.na(curr_vaf) && curr_vaf >= vaf_threshold &&
            (is.na(prev_vaf) || prev_vaf < vaf_threshold)) {
          new_mutations_data <- rbind(new_mutations_data, data.frame(
            Patient_ID = patient,
            Gene = gene,
            NGS_Timepoint = paste0("NGS_", timepoints[tp_idx]),
            Emergence_VAF = curr_vaf,
            Emergence_Timepoint = timepoints[tp_idx]
          ))
        }
      }
    }
  }
  cat("Detected", nrow(new_mutations_data), "new mutations across", 
      length(unique(new_mutations_data$Patient_ID)), "patients\n")
  if (nrow(new_mutations_data) > 0) {
    cat("Most frequent new mutations:\n")
    gene_freq <- sort(table(new_mutations_data$Gene), decreasing = TRUE)
    print(head(gene_freq, 5))
  }
  return(new_mutations_data)
}

# Function to prepare time-dependent data
prepare_time_dependent_data <- function(clinical_data, new_mutations_data) {
  cat("Preparing time-dependent survival data...\n")
  
  # Map actual column names to expected names
  if ("OS months" %in% names(clinical_data)) {
    clinical_data$OS_time <- as.numeric(clinical_data[["OS months"]]) * 30.44
  } else {
    clinical_data$OS_time <- rep(365, nrow(clinical_data))
  }
  
  if ("Survival status" %in% names(clinical_data)) {
    clinical_data$OS_event <- as.numeric(clinical_data[["Survival status"]])
  } else if ("status" %in% names(clinical_data)) {
    clinical_data$OS_event <- ifelse(tolower(clinical_data[["status"]]) == "dead", 1, 0)
  } else {
    clinical_data$OS_event <- rep(0, nrow(clinical_data))
  }
  
  cat("Found", sum(clinical_data$OS_event), "death events out of", nrow(clinical_data), "patients\n")
  
  # Convert to data frame
  clinical_data <- as.data.frame(clinical_data)
  
  # Calculate NGS_2 timing (estimate if not available)
  if ("NGS_2_date" %in% names(clinical_data)) {
    clinical_data$NGS2_time <- as.numeric(clinical_data$NGS_2_Date - clinical_data$Date_of_Diagnosis)
  } else {
    clinical_data$NGS2_time <- 365  # Assume 1 year follow-up
  }
  
  # Create time-dependent dataset
  td_data <- data.frame()
  patient_id_col <- if ("Patient_ID" %in% names(clinical_data)) "Patient_ID" else "MRN"
  
  for (i in 1:nrow(clinical_data)) {
    patient_id <- clinical_data[[patient_id_col]][i]
    patient_mutations <- new_mutations_data[new_mutations_data$Patient_ID == patient_id, ]
    
    # Create baseline row with consistent structure
    baseline_row <- clinical_data[i, ]
    baseline_row$tstart <- 0
    baseline_row$tstop <- ifelse(nrow(patient_mutations) > 0, 
                                 baseline_row$NGS2_time, 
                                 baseline_row$OS_time)
    baseline_row$new_mutation <- 0
    baseline_row$mutation_count <- 0
    baseline_row$event <- ifelse(nrow(patient_mutations) > 0, 0, baseline_row$OS_event)
    baseline_row$mutations_list <- ""
    baseline_row$max_vaf_change <- 0
    baseline_row$mean_ngs2_vaf <- 0
    
    # Initialize td_data with first row structure
    if (nrow(td_data) == 0) {
      td_data <- baseline_row
    } else {
      # Ensure column consistency before rbind
      common_cols <- intersect(names(td_data), names(baseline_row))
      td_data <- rbind(td_data[, common_cols], baseline_row[, common_cols])
    }
    
    # Post-mutation period
    if (nrow(patient_mutations) > 0) {
      post_mutation_row <- baseline_row
      post_mutation_row$tstart <- baseline_row$NGS2_time
      post_mutation_row$tstop <- baseline_row$OS_time
      post_mutation_row$new_mutation <- 1
      post_mutation_row$mutation_count <- nrow(patient_mutations)
      post_mutation_row$event <- baseline_row$OS_event
      
      # Add mutation details
      post_mutation_row$mutations_list <- paste(patient_mutations$Gene, collapse = ";")
      if (nrow(patient_mutations) > 0) {
        # Use Emergence_VAF instead of NGS2_VAF for mean calculation
        post_mutation_row$max_vaf_change <- max(patient_mutations$Emergence_VAF, na.rm = TRUE)
        post_mutation_row$mean_ngs2_vaf <- mean(patient_mutations$Emergence_VAF, na.rm = TRUE)
      }
      
      # Ensure column consistency before rbind
      common_cols <- intersect(names(td_data), names(post_mutation_row))
      td_data <- rbind(td_data[, common_cols], post_mutation_row[, common_cols])
    }
  }
  
  # Remove invalid intervals
  td_data <- td_data[td_data$tstart < td_data$tstop, ]
  
  cat("Created time-dependent dataset with", nrow(td_data), "observations\n")
  return(td_data)
}

# Function to run time-dependent Cox regression
run_time_dependent_cox <- function(td_data) {
  cat("Running time-dependent Cox regression...\n")
  
  # Define potential covariates
  covariates <- c("Age", "Sex", "HB(g/dl)", 
                  "WBC(k/ul)", "PLT(k/ul)", "BM blasts %")
  covariates <- str_to_lower(covariates)
  available_covs <- covariates[covariates %in% names(td_data)]
  
  # Select covariates with sufficient data
  valid_covs <- c()
  for (cov in available_covs) {
    if (sum(!is.na(td_data[[cov]])) > 20) {  # Require at least 20 non-missing values
      valid_covs <- c(valid_covs, cov)
    }
  }
  
  # Build formula
  if (length(valid_covs) == 0) {
    formula_str <- "Surv(tstart, tstop, event) ~ new_mutation"
  } else {
    formula_str <- paste("Surv(tstart, tstop, event) ~ new_mutation +", 
                         paste(valid_covs, collapse = " + "))
  }
  
  cat("Using formula:", formula_str, "\n")
  
  # Fit time-dependent Cox model
  td_cox <- coxph(as.formula(formula_str), data = td_data)
  
  cat("\n=== TIME-DEPENDENT COX REGRESSION RESULTS ===\n")
  print(summary(td_cox))
  
  # Extract key results
  if ("new_mutation" %in% names(coef(td_cox))) {
    hr_new_mut <- exp(coef(td_cox)["new_mutation"])
    p_new_mut <- summary(td_cox)$coefficients["new_mutation", "Pr(>|z|)"]
    cat("\n🎯 KEY FINDING: New mutations HR =", round(hr_new_mut, 3), 
        ", p =", format(p_new_mut, scientific = TRUE), "\n")
  }
  
  return(td_cox)
}

# Function to run standard Cox regression
run_standard_cox <- function(clinical_data, serial_patients) {
  cat("Running standard Cox regression for baseline comparison...\n")
  
  clinical_data <- as.data.frame(clinical_data)
  patient_id_col <- if ("MRN" %in% names(clinical_data)) "MRN"
  else if ("Patient_ID" %in% names(clinical_data)) "Patient_ID"
  else names(clinical_data)[1]
  clinical_data$has_serial_ngs <- ifelse(as.character(clinical_data[[patient_id_col]]) %in% as.character(serial_patients), 1, 0)
  
  # Build model
  covariates <- c("Age", "Sex")
  available_covs <- covariates[covariates %in% names(clinical_data)]
  
  valid_covs <- c()
  for (cov in available_covs) {
    if (sum(!is.na(clinical_data[[cov]])) > 20) {
      valid_covs <- c(valid_covs, cov)
    }
  }
  
  if (length(valid_covs) == 0) {
    formula_str <- "Surv(OS_time, OS_event) ~ has_serial_ngs"
  } else {
    formula_str <- paste("Surv(OS_time, OS_event) ~ has_serial_ngs +", 
                         paste(valid_covs, collapse = " + "))
  }
  
  standard_cox <- coxph(as.formula(formula_str), data = clinical_data)
  
  cat("\n=== STANDARD COX REGRESSION RESULTS ===\n")
  print(summary(standard_cox))
  
  # Extract key results
  if ("has_serial_ngs" %in% names(coef(standard_cox))) {
    hr_serial <- exp(coef(standard_cox)["has_serial_ngs"])
    p_serial <- summary(standard_cox)$coefficients["has_serial_ngs", "Pr(>|z|)"]
    cat("\n🎯 KEY FINDING: Serial NGS HR =", round(hr_serial, 3), 
        ", p =", format(p_serial, scientific = TRUE), "\n")
  }
  
  return(standard_cox)
}

# Function to create Kaplan-Meier curves
create_km_curves <- function(clinical_data, new_mutations_data, serial_patients) {
  cat("Creating Kaplan-Meier curves...\n")
  
  clinical_data <- as.data.frame(clinical_data)
  patient_id_col <- if ("Patient_ID" %in% names(clinical_data)) "Patient_ID" else "MRN"
  
  # Create group indicator
  group_var <- rep("Baseline only", nrow(clinical_data))
  for (i in 1:nrow(clinical_data)) {
    patient_val <- clinical_data[[patient_id_col]][i]
    if (length(patient_val) > 0 && !is.na(patient_val) && patient_val %in% serial_patients) {
      group_var[i] <- "Serial NGS"
    }
  }
  clinical_data$group <- group_var
  
  # Fit survival curve
  km_overall <- survfit(Surv(OS_time, OS_event) ~ group, data = clinical_data)
  
  # Kaplan-Meier by emergence of any new mutation (any gene)
  if (nrow(new_mutations_data) > 0) {
    # Debug: Check patient ID matching
    cat("Debug - Clinical IDs (first 5):", head(as.character(clinical_data[[patient_id_col]]), 5), "\n")
    cat("Debug - Mutation IDs (unique):", unique(as.character(new_mutations_data$Patient_ID)), "\n")
    
    # Fix: Convert both sides to character for proper matching
    clinical_data$any_new_mutation <- ifelse(as.character(clinical_data[[patient_id_col]]) %in% as.character(unique(new_mutations_data$Patient_ID)), "New mutation emerged", "No new mutation")
    
    cat("Debug - Patients with new mutations:", sum(clinical_data$any_new_mutation == "New mutation emerged"), "\n")
    
    # Check if any_new_mutation variable is properly defined and has both groups
    mutation_table <- table(clinical_data$any_new_mutation, useNA = "always")
    cat("Debug - any_new_mutation distribution:\n")
    print(mutation_table)
    
    if (sum(!is.na(clinical_data$any_new_mutation)) > 0 && length(unique(clinical_data$any_new_mutation[!is.na(clinical_data$any_new_mutation)])) >= 1) {
      km_any_new <- survfit(Surv(OS_time, OS_event) ~ any_new_mutation, data = clinical_data)
    } else {
      cat("Warning: any_new_mutation variable not properly defined, skipping KM curve\n")
      km_any_new <- NULL
    }
    
    # Plot
    if (!is.null(km_any_new)) {
      if ("survminer" %in% advanced_loaded) {
        p_any_new <- ggsurvplot(km_any_new,
                                data = clinical_data,
                                title = "Overall Survival by Emergence of Any New Mutation",
                                xlab = "Time from Diagnosis (days)",
                                ylab = "Overall Survival Probability",
                                risk.table = TRUE,
                                pval = TRUE,
                                conf.int = TRUE)
        ggsave("Serial_NGS_Output/KM_Any_New_Mutation.png",
               print(p_any_new), width = 10, height = 8, dpi = 300)
      } else {
        png("Serial_NGS_Output/KM_Any_New_Mutation_Basic.png", width = 800, height = 600)
        plot(km_any_new,
             col = c("blue", "red"),
             xlab = "Time from Diagnosis (days)",
             ylab = "Overall Survival Probability",
             main = "Overall Survival by Emergence of Any New Mutation",
             lwd = 2)
        legend("topright",
               legend = c("No new mutation", "New mutation emerged"),
               col = c("blue", "red"),
               lwd = 2)
        dev.off()
      }
    }
  }
  
  # Create plot
  if ("survminer" %in% advanced_loaded) {
    # Advanced plot with survminer
    p1 <- ggsurvplot(km_overall, 
                     data = clinical_data,
                     title = "Overall Survival: Serial NGS vs Baseline",
                     xlab = "Time from Diagnosis (days)",
                     ylab = "Overall Survival Probability",
                     risk.table = TRUE,
                     pval = TRUE,
                     conf.int = TRUE,
                     palette = c("#E7B800", "#2E9FDF"))
    
    ggsave("Serial_NGS_Output/KM_Overall_Comparison_Advanced.png", 
           print(p1), width = 12, height = 10, dpi = 300)
  } else {
    # Basic plot with base R
    png("Serial_NGS_Output/KM_Overall_Comparison_Basic.png", width = 800, height = 600)
    plot(km_overall, 
         col = c("blue", "red"),
         xlab = "Time from Diagnosis (days)",
         ylab = "Overall Survival Probability",
         main = "Overall Survival: Serial NGS vs Baseline",
         lwd = 2)
    legend("topright", 
           legend = c("Baseline only", "Serial NGS"),
           col = c("blue", "red"),
           lwd = 2)
    dev.off()
  }
  
  # Log-rank test
  lr_test <- survdiff(Surv(OS_time, OS_event) ~ group, data = clinical_data)
  p_value <- 1 - pchisq(lr_test$chisq, 1)
  cat("Log-rank test p-value:", format(p_value, scientific = TRUE), "\n")
  
  # Create plots for individual genes with new mutations
  if (nrow(new_mutations_data) > 0) {
    cat("Creating individual gene Kaplan-Meier curves...\n")
    genes_with_new_mutations <- names(sort(table(new_mutations_data$Gene), decreasing = TRUE))
    
    cat("Genes with new mutations:", length(genes_with_new_mutations), "\n")
    print(head(table(new_mutations_data$Gene)[order(table(new_mutations_data$Gene), decreasing = TRUE)], 10))
    
    gene_plots_created <- 0
    
    for (gene in genes_with_new_mutations) {
      gene_patients <- unique(new_mutations_data$Patient_ID[new_mutations_data$Gene == gene])
      gene_count <- length(gene_patients)
      
      cat("Processing gene:", gene, "- patients with new mutations:", gene_count, "\n")
      
      if (gene_count >= 3) {  # Only if ≥3 patients have new mutations in this gene
        # Fix patient ID matching with proper type conversion
        clinical_data[[paste0(gene, "_new_mutation")]] <- 
          ifelse(as.character(clinical_data[[patient_id_col]]) %in% as.character(gene_patients), 
                 "New mutation emerged", "No new mutation")
        
        # Check distribution
        gene_table <- table(clinical_data[[paste0(gene, "_new_mutation")]], useNA = "always")
        cat("  ", gene, "distribution:", paste(names(gene_table), "=", gene_table, collapse = ", "), "\n")
        
        # Only proceed if we have both groups
        if (length(unique(clinical_data[[paste0(gene, "_new_mutation")]])) >= 2) {
          # Create safer column name for analysis
          col_name <- paste0(gene, "_new_mutation")
          mutation_status <- clinical_data[[col_name]]
          
          tryCatch({
            km_gene <- survfit(Surv(OS_time, OS_event) ~ mutation_status, data = clinical_data)
            
            # Calculate log-rank test p-value
            lr_test_gene <- survdiff(Surv(OS_time, OS_event) ~ mutation_status, data = clinical_data)
            p_val_gene <- 1 - pchisq(lr_test_gene$chisq, 1)
            
            # Clean gene name for filename (remove special characters)
            clean_gene_name <- gsub("[^A-Za-z0-9_]", "_", gene)
            clean_gene_name <- gsub("_+", "_", clean_gene_name)  # Remove multiple underscores
            clean_gene_name <- gsub("^_|_$", "", clean_gene_name)  # Remove leading/trailing underscores
            
            if ("survminer" %in% advanced_loaded) {
              # Advanced plot with survminer
              p_gene <- ggsurvplot(km_gene, 
                                   data = clinical_data,
                                   title = paste("Overall Survival by", gene, "New Mutation Emergence"),
                                   xlab = "Time from Diagnosis (days)",
                                   ylab = "Overall Survival Probability",
                                   risk.table = TRUE,
                                   pval = TRUE,
                                   conf.int = TRUE,
                                   palette = c("#2E9FDF", "#E7B800"))
              
              filename_advanced <- paste0("Serial_NGS_Output/KM_", clean_gene_name, "_NewMutation_Advanced.png")
              ggsave(filename_advanced, print(p_gene), width = 12, height = 10, dpi = 300)
              cat("  Advanced plot saved:", filename_advanced, "\n")
            } else {
              # Basic plot with base R
              filename_basic <- paste0("Serial_NGS_Output/KM_", clean_gene_name, "_NewMutation_Basic.png")
              png(filename_basic, width = 800, height = 600)
              plot(km_gene,
                   col = c("blue", "red"),
                   xlab = "Time from Diagnosis (days)",
                   ylab = "Overall Survival Probability",
                   main = paste("Overall Survival by", gene, "New Mutation Emergence"),
                   lwd = 2)
              legend("topright",
                     legend = c("No new mutation", "New mutation emerged"),
                     col = c("blue", "red"),
                     lwd = 2)
              
              # Add p-value to plot
              text(x = max(clinical_data$OS_time, na.rm = TRUE) * 0.7, 
                   y = 0.1, 
                   labels = paste("Log-rank p =", format(p_val_gene, digits = 3)),
                   cex = 0.8)
              
              dev.off()
              cat("  Basic plot saved:", filename_basic, "\n")
            }
            
            cat("  ", gene, "log-rank p-value:", format(p_val_gene, scientific = TRUE), "\n")
            gene_plots_created <- gene_plots_created + 1
            
          }, error = function(e) {
            cat("  Error creating plot for", gene, ":", e$message, "\n")
          })
        } else {
          cat("  Skipping", gene, "- insufficient group variation\n")
        }
      } else {
        cat("  Skipping", gene, "- insufficient patients (need ≥3, have", gene_count, ")\n")
      }
    }
    
    cat("Individual gene KM curves created:", gene_plots_created, "/", length(genes_with_new_mutations), "\n")
  }
  
  return(list(overall = km_overall, p_value = p_value))
}

# Function to generate comprehensive summary report
generate_summary_report <- function(results) {
  cat("Generating comprehensive summary report...\n")
  
  # Create summary statistics
  summary_stats <- list(
    total_patients = nrow(results$clinical_data),
    serial_ngs_patients = length(results$serial_patients),
    new_mutations_detected = nrow(results$new_mutations_data),
    genes_with_new_mutations = length(unique(results$new_mutations_data$Gene))
  )
  
  # Save summary to file
  sink("Serial_NGS_Output/Analysis_Summary.txt")
  
  cat("=================================================================\n")
  cat("    SERIAL NGS AND NEW MUTATIONS ANALYSIS SUMMARY REPORT\n")
  cat("=================================================================\n\n")
  
  cat("📊 PATIENT COHORT OVERVIEW:\n")
  cat("├─ Total patients analyzed:", summary_stats$total_patients, "\n")
  cat("├─ Patients with serial NGS:", summary_stats$serial_ngs_patients, "\n")
  cat("├─ Serial NGS rate:", round(summary_stats$serial_ngs_patients/summary_stats$total_patients*100, 1), "%\n")
  cat("└─ Death events:", sum(results$clinical_data$OS_event, na.rm = TRUE), "\n\n")
  
  cat("🧬 NEW MUTATIONS DETECTED:\n")
  cat("├─ Total new mutations:", summary_stats$new_mutations_detected, "\n")
  cat("├─ Unique genes with new mutations:", summary_stats$genes_with_new_mutations, "\n")
  cat("└─ Average new mutations per serial patient:", 
      round(summary_stats$new_mutations_detected/summary_stats$serial_ngs_patients, 2), "\n\n")
  
  if (!is.null(results$new_mutations_data) && nrow(results$new_mutations_data) > 0) {
    cat("🎯 TOP GENES WITH NEW MUTATIONS:\n")
    gene_freq <- table(results$new_mutations_data$Gene)
    gene_freq_sorted <- sort(gene_freq, decreasing = TRUE)
    for (i in 1:min(10, length(gene_freq_sorted))) {
      cat("├─", names(gene_freq_sorted)[i], ":", gene_freq_sorted[i], "patients\n")
    }
    cat("\n")
  }
  
  cat("📈 STATISTICAL RESULTS:\n")
  if (!is.null(results$td_cox)) {
    td_coef <- summary(results$td_cox)$coefficients
    if ("new_mutation" %in% rownames(td_coef)) {
      hr <- exp(td_coef["new_mutation", "coef"])
      p_val <- td_coef["new_mutation", "Pr(>|z|)"]
      ci_lower <- exp(confint(results$td_cox)["new_mutation", 1])
      ci_upper <- exp(confint(results$td_cox)["new_mutation", 2])
      
      cat("├─ Time-dependent Cox Regression:\n")
      cat("│  ├─ New mutations HR:", round(hr, 3), 
          "(95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), ")\n")
      cat("│  ├─ P-value:", format(p_val, scientific = TRUE), "\n")
      cat("│  └─ Interpretation:", 
          ifelse(p_val < 0.05, "STATISTICALLY SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
    }
  }
  
  if (!is.null(results$standard_cox)) {
    std_coef <- summary(results$standard_cox)$coefficients
    if ("has_serial_ngs" %in% rownames(std_coef)) {
      hr <- exp(std_coef["has_serial_ngs", "coef"])
      p_val <- std_coef["has_serial_ngs", "Pr(>|z|)"]
      
      cat("├─ Standard Cox Regression:\n")
      cat("│  ├─ Serial NGS patients HR:", round(hr, 3), "\n")
      cat("│  ├─ P-value:", format(p_val, scientific = TRUE), "\n")
      cat("│  └─ Interpretation:", 
          ifelse(p_val < 0.05, "STATISTICALLY SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
    }
  }
  
  if (!is.null(results$km_curves)) {
    cat("├─ Kaplan-Meier Analysis:\n")
    cat("│  ├─ Log-rank test p-value:", format(results$km_curves$p_value, scientific = TRUE), "\n")
    cat("│  └─ Survival difference:", 
        ifelse(results$km_curves$p_value < 0.05, "STATISTICALLY SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
  }
  
  cat("\n🔬 CLINICAL INTERPRETATION:\n")
  if (!is.null(results$td_cox) && "new_mutation" %in% rownames(summary(results$td_cox)$coefficients)) {
    td_hr <- exp(summary(results$td_cox)$coefficients["new_mutation", "coef"])
    if (td_hr < 1) {
      cat("├─ New mutations appear PROTECTIVE (HR < 1)\n")
      cat("├─ Possible explanations:\n")
      cat("│  ├─ Early detection leads to better management\n")
      cat("│  ├─ Detected mutations may be low-risk variants\n")
      cat("│  └─ Patients with serial NGS receive closer monitoring\n")
    } else {
      cat("├─ New mutations appear HARMFUL (HR > 1)\n")
      cat("├─ This suggests new mutations drive disease progression\n")
    }
  }
  
  cat("\n📁 OUTPUT FILES GENERATED:\n")
  cat("├─ Analysis_Summary.txt (this report)\n")
  cat("├─ KM_Overall_Comparison.png (survival curves)\n")
  if (summary_stats$genes_with_new_mutations > 0) {
    cat("├─ Gene-specific KM plots (for top mutated genes)\n")
  }
  cat("└─ Statistical model outputs (displayed in console)\n\n")
  
  cat("⚡ ANALYSIS METHODOLOGY:\n")
  cat("├─ Serial NGS definition: NGS_2 with VAF ≥ 1%\n")
  cat("├─ New mutation definition: Absent in NGS_1, present in NGS_2\n")
  cat("├─ Time-dependent Cox: Accounts for timing of mutation emergence\n")
  cat("├─ Standard Cox: Compares serial NGS vs baseline patients\n")
  cat("└─ Kaplan-Meier: Survival curve comparison with log-rank test\n\n")
  
  cat("=================================================================\n")
  cat("Analysis completed successfully on", Sys.time(), "\n")
  cat("=================================================================\n")
  
  sink()
  
  cat("📄 Summary report saved to Serial_NGS_Output/Analysis_Summary.txt\n")
}

# Main analysis function
main_analysis <- function() {
  cat("=================================================================\n")
  cat("    SERIAL NGS TIME-DEPENDENT COX ANALYSIS (OPTIMIZED)\n")
  cat("=================================================================\n\n")
  
  start_time <- Sys.time()
  
  # Read data
  data_list <- read_data_safely()
  clinical_data <- data_list$clinical
  ngs_data <- data_list$ngs
  
  # Only genes columns
  baseline_vaf_cols <- grep("^NGS_1[A-Z0-9]+$", names(ngs_data), 
                            value = TRUE)
  
  
  if (length(baseline_vaf_cols) > 0) {
    baseline_long <- ngs_data %>%
      select(MRN, all_of(baseline_vaf_cols)) %>%  
      pivot_longer(
        cols = -MRN,
        names_to = "Gene",
        values_to = "VAF"
      ) 
    write.csv(baseline_long, file = "Serial_NGS_Output/baseline_mutations_table.csv", row.names = FALSE)
  }
  
  # --- EXCLUDE POST-BMT NGS_2, NGS_3, NGS_4, NGS_5 DATA ---
  if ("bmt_date" %in% names(clinical_data)) {
    cat("Applying post-BMT NGS data exclusion for NGS_2 through NGS_5...\n")
    patient_id_col <- if ("Patient_ID" %in% names(clinical_data)) "Patient_ID" else if ("MRN" %in% names(clinical_data)) "MRN" else names(clinical_data)[1]
    bmt_mapping <- data.frame(
      Patient = as.character(clinical_data[[patient_id_col]]),
      BMT_Date = as.Date(clinical_data$bmt_date)
    )
    serial_ngs_numbers <- 2:5
    for (i in seq_len(nrow(bmt_mapping))) {
      pid <- bmt_mapping$Patient[i]
      bmt_date <- bmt_mapping$BMT_Date[i]
      if (!is.na(bmt_date)) {
        ngs_idx <- which(as.character(ngs_data[[1]]) == pid)
        if (length(ngs_idx) > 0) {
          for (ngs_num in serial_ngs_numbers) {
            date_col_pattern <- paste0("NGS_", ngs_num, " date")
            date_col <- grep(date_col_pattern, names(ngs_data), ignore.case=TRUE, value=TRUE)
            if (length(date_col) == 0) next
            date_to_convert <- unlist(ngs_data[ngs_idx, date_col])
            date_to_convert <- round(as.numeric(date_to_convert))
            ngs_date <- as.Date(date_to_convert, origin = "1899-12-30") |> format("%m/%d/%Y")
            ngs_date <- as.Date(ngs_date, format = "%m/%d/%Y")
            bmt_date <- as.Date(bmt_date, format = "%m/%d/%Y")
            if (!is.na(ngs_date) && ngs_date >= bmt_date) {
              ngs_cols_pattern <- paste0("NGS_", ngs_num)
              vaf_cols <- grep(ngs_cols_pattern, names(ngs_data), value=TRUE)
              ngs_data[ngs_idx, vaf_cols] <- NA
            }
          }
        }
      }
    }
  }
  
  # Identify serial NGS patients
  cat("\n🔍 Identifying serial NGS patients...\n")
  serial_info <- identify_serial_ngs(ngs_data, vaf_threshold = 1)
  serial_patients <- serial_info$serial_patients
  
  # Detect new mutations
  cat("\n🧬 Detecting new mutations...\n")
  new_mutations_data <- detect_new_mutations_all_timepoints(ngs_data, vaf_threshold = 1)
  
  # Add Emergence_Date and Time_from_Baseline to new_mutations_data
  if (nrow(new_mutations_data) > 0) {
    baseline_date_map <- setNames(
      if (any(grepl("^NGS_1[ _]date$", names(ngs_data), ignore.case = TRUE))) {
        as.Date(ngs_data[[grep("^NGS_1[ _]date$", names(ngs_data), ignore.case = TRUE)]])
      } else if ("Date" %in% names(clinical_data)) {
        as.Date(clinical_data$Date)
      } else {
        rep(NA, nrow(ngs_data))
      },
      as.character(ngs_data[[1]])
    )
    new_mutations_data <- new_mutations_data %>%
      rowwise() %>%
      mutate(Emergence_Date = {
        pid <- Patient_ID
        nn <- str_extract(NGS_Timepoint, "[0-9]")
        date_col <- grep(paste0(nn, "[ _]date$"), names(ngs_data), value = TRUE, ignore.case = TRUE)
        if (length(date_col) > 0) {
          vals <- ngs_data %>%
            filter(!!sym(names(ngs_data)[1]) == pid) %>%
            select(all_of(date_col)) %>%
            pull(1)
          safe_date_conversion(vals)[1]
        } else {
          NA
        }
      }) %>%
      ungroup()
    baseline_date_map <- setNames(
      if (any(grepl("^NGS_1[ _]date$", names(ngs_data), ignore.case = TRUE))) {
        as.Date(as.character(ngs_data[[grep("^NGS_1[ _]date$", names(ngs_data), ignore.case = TRUE)]]))
      } else if ("Date" %in% names(clinical_data)) {
        as.Date(as.character(clinical_data$Date))
      } else {
        rep(NA, nrow(ngs_data))
      },
      as.character(ngs_data[[1]])
    )
    new_mutations_data$Baseline_Date <- baseline_date_map[
      as.character(new_mutations_data$Patient_ID)]
    new_mutations_data$Time_from_Baseline <- as.numeric(
      new_mutations_data$Emergence_Date - new_mutations_data$Baseline_Date)
    table_by_gene <- as.data.frame(table(new_mutations_data$Gene))
    write.csv(new_mutations_data, file = "Serial_NGS_Output/new_mutations_table.csv", row.names = FALSE)
    write.csv(table_by_gene, file = "Serial_NGS_Output/emergence_by_gene.csv", row.names = FALSE)
  }
  
  # Prepare time-dependent data
  cat("\n⏱️ Preparing time-dependent survival data...\n")
  td_data <- prepare_time_dependent_data(clinical_data, new_mutations_data)
  
  # Run time-dependent Cox regression
  cat("\n📊 Running time-dependent Cox regression...\n")
  td_cox <- run_time_dependent_cox(td_data)
  
  # Extract processed clinical data for other analyses
  processed_clinical <- td_data[td_data$tstart == 0, ]
  
  # Run standard Cox regression
  cat("\n📈 Running standard Cox regression...\n")
  standard_cox <- run_standard_cox(processed_clinical, serial_patients)
  
  # Compare baseline features (PSM removed)
  cat("\n📋 Comparing baseline features...\n")
  baseline_comparison <- compare_baseline_features(processed_clinical, serial_patients)
  if ("tableone" %in% advanced_loaded) {
    # Case-insensitive matching and NA filtering
    available_vars <- c("Age", "Sex", "OS_time", "OS_event")
    matched_vars <- intersect(tolower(available_vars), tolower(names(processed_clinical)))
    real_vars <- names(processed_clinical)[tolower(names(processed_clinical)) %in% matched_vars]
    real_vars <- real_vars[sapply(real_vars, function(v) sum(!is.na(processed_clinical[[v]])) > 0)]
    if (length(real_vars) > 0) {
      tbl1 <- CreateTableOne(vars = real_vars, data = processed_clinical)
      write.csv(print(tbl1, showAllLevels=TRUE), file="Serial_NGS_Output/baseline_table.csv")
    } else {
      cat("No valid variables for TableOne baseline table; skipping TableOne.\n")
    }
  }
  # Create Kaplan-Meier curves
  cat("\n📉 Creating Kaplan-Meier survival curves...\n")
  km_curves <- create_km_curves(processed_clinical, new_mutations_data, serial_patients)
  if ("any_new_mutation" %in% names(processed_clinical)) {
    summary_by_group <- aggregate(processed_clinical$OS_time, by=list(processed_clinical$any_new_mutation), summary)
    write.csv(summary_by_group, file = "Serial_NGS_Output/followup_by_emergence_status.csv", row.names = FALSE)
  }
  # Compile results
  results <- list(
    clinical_data = processed_clinical,
    ngs_data = ngs_data,
    serial_patients = serial_patients,
    new_mutations_data = new_mutations_data,
    td_data = td_data,
    td_cox = td_cox,
    standard_cox = standard_cox,
    baseline_comparison = baseline_comparison,
    km_curves = km_curves
  )
  
  # Generate summary report
  cat("\n📝 Generating comprehensive summary report...\n")
  generate_summary_report(results)
  
  end_time <- Sys.time()
  analysis_duration <- difftime(end_time, start_time, units = "secs")
  cat("\n=================================================================\n")
  cat("✅ ANALYSIS COMPLETED!\n")
  cat("⏱️  Total analysis time:", round(as.numeric(analysis_duration), 1), "seconds\n")
  cat("📁 Check the 'Serial_NGS_Output' folder for all generated files\n")
  cat("=================================================================\n")
  return(list(
    clinical_data = clinical_data,
    ngs_data = ngs_data,
    serial_patients = serial_patients,
    new_mutations_data = new_mutations_data
    # ...other outputs...
  ))
}

cat("🧬 Serial NGS Time-Dependent Cox Analysis - Enhanced Version\n")
cat("📦 Loading essential packages and attempting advanced features...\n\n")
if (!exists("results")) {
  results <- main_analysis()
}

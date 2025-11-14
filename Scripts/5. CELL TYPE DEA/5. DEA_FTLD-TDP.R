library(edgeR, lib.loc = "/home/jaumatell/R/x86_64-pc-linux-gnu-library/4.4/edgeR_4")
library("dplyr")
library("GO.db")
library("xlsx")
################################################################################
# Paths and parameters
CSV_DECONVOLDED_PATH <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE/"
CASE_LEGEND <- "/media/jaumatell/datos/URI/BayesPrism/FTD/DATASETS/FTD_QUIM/decoder_DeSeq2_FTD_FINAL.xlsx"
OUTPUT_DIRECTORY <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/FTLD/TDP/DEG"
REFERENCE_GROUP <- "Healthy"
################################################################################

csv_files <- list.files(CSV_DECONVOLDED_PATH, pattern = "\\.csv$", full.names = TRUE)

# Read case legend and prepare group IDs
case_legend <- read.xlsx(CASE_LEGEND,sheetIndex = 1)
case_legend$group_ID <- factor(case_legend$group.ID)
case_legend$group.ID <- NULL

# Ensure levels are consistent and drop unused levels
conditions <- levels(as.factor(case_legend$group_ID))
conditions <- "TDP"

# Create output directory if it doesn't exist
if (!file.exists(OUTPUT_DIRECTORY)) {
  dir.create(OUTPUT_DIRECTORY, recursive = TRUE)
}

# Initialize summary matrix
matrix_data <- matrix(1:length(conditions)*length(csv_files), nrow = length(conditions), ncol = length(csv_files))
df <- data.frame(matrix_data, row.names = conditions)
colnames(df) <- sapply(csv_files, function(file) tools::file_path_sans_ext(basename(file)))

# Process each CSV file
for (file in csv_files) {
  filename <- tools::file_path_sans_ext(basename(file))
  print(filename)
  
  data <- read.csv(file, row.names = 1)
  rownames(data) <-gsub("X", "", rownames(data))
  celltype_folder <- file.path(OUTPUT_DIRECTORY, filename)
  if (!file.exists(celltype_folder)) {
    dir.create(celltype_folder, recursive = TRUE)
  }
  
  for (condition in conditions) {
    condition_folder <- file.path(celltype_folder, condition)
    if (!file.exists(condition_folder)) {
      dir.create(condition_folder, recursive = TRUE)
    }
    
    tryCatch({
      results_table <- NULL
      res <- NULL
      
      comparison_data <- case_legend[case_legend$group_ID %in% c(REFERENCE_GROUP, condition), ]
      comparison_data$group_ID <- relevel(factor(comparison_data$group_ID, levels = c(REFERENCE_GROUP, condition)), ref = condition)
      comparison_data <- comparison_data[order(comparison_data$group_ID), ]
      
      common_ids <- intersect(comparison_data$sample.ID, rownames(data))
      data_subset <- data[rownames(data) %in% common_ids, ]
      data_ordered <- data_subset[match(comparison_data$sample.ID, rownames(data_subset)), ]
      data_ordered[is.na(data_ordered)] <- 10e-8
      #data_ordered <- data_ordered[row_sums > 0, ]
      
      group <- factor(comparison_data$group_ID)
      y <- DGEList(counts = t(data_ordered), group = group)
      
      # Plot average log CPM
      AveLogCPM <- aveLogCPM(y)
      y <- normLibSizes(y)
      #y <- y[row_sums > 0, ]
      
      plot_filename <- file.path(condition_folder, "MD_plot_h.png")
      png(file = plot_filename)
      plotMD(y, column = 1)
      abline(h = 0, col = "red", lty = 2, lwd = 2)
      dev.off()
      
      # Create design matrix
      design <- model.matrix(~ 0 + group)
      colnames(design) <- levels(group)
      
      # Estimate dispersion
      y <- estimateDisp(y, design)
      #y <- EdgeR::fillterByExpr(y, group=group)
      fit <- glmQLFit(y, design)
      
      #Differential expression analysis
      contrast <- makeContrasts(contrasts = paste0(make.names(condition), " - ", make.names(REFERENCE_GROUP)), levels = design)
      res <- glmLRT(fit, contrast = contrast)
      is.de <- decideTests(res, adjust.method = "fdr", p.value = 0.05, lfc = 0)
      results_table <- res$table
      results_table$adj_pval <- p.adjust(results_table$PValue, method = "fdr")
      results_table <- cbind(results_table, data.frame(is.de))
      
      # Save results
      write.csv(res$table, file = file.path(condition_folder, "res_h.csv"))
      write.csv(results_table, file = file.path(condition_folder, "results_adj_h.csv"))
      
      # Histogram
      plot_filename <- file.path(condition_folder, "histogram_plot_h.png")
      png(file = plot_filename)
      hist(AveLogCPM)
      dev.off()
      
      # BCV plot
      plot_filename <- file.path(condition_folder, "BCV_plot_h.png")
      png(file = plot_filename)
      plotBCV(y)
      dev.off()
      
      # "Volcano like" plot
      plot_filename <- file.path(condition_folder, "MD_res_plot_h.png")
      png(file = plot_filename)
      plotMD(res, status = is.de)
      dev.off()
      
      # Heatmap clustering
      logCPM <- cpm(y, prior.count = 2, log = TRUE)
      colnames(logCPM) <- paste(y$samples$group, 1:2, sep = "-")
      tr <- glmTreat(fit, contrast = contrast, lfc = log2(1.5))
      o <- order(tr$table$PValue)
      logCPM <- logCPM[o[1:30], ]
      
      plot_filename <- file.path(condition_folder, "Heatmap_h.png")
      png(file = plot_filename)
      coolmap(logCPM, margins = c(7, 7), lhei = c(1, 6), lwid = c(1, 3))
      dev.off()
      
      # Update summary DataFrame
      df[condition, filename] <- nrow(results_table[results_table$adj_pval < 0.05, ])
      
    }, error = function(e) {
      print(e)
      df[condition, filename] <- 0
    })
  }
}

# Save summary DataFrame
write.csv(df, file = file.path(OUTPUT_DIRECTORY, "summary_h.csv"))

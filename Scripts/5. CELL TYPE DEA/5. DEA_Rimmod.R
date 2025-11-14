library("edgeR")
library("dplyr")
library("GO.db")
library("xlsx")
################################################################################
# Paths and parameters
CSV_DECONVOLDED_PATH <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/RIMOD/CELL STATE ORIGINAL/"
CASE_LEGEND <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/RIMOD_BULK/DATA/rimod_ftd_dataset_table_v3.txt"
OUTPUT_DIRECTORY <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/RIMOD/CS_SV"
SV_COVARIATE <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/RIMOD_BULK/METADATA/SVA.csv"

REFERENCE_GROUP <- "Healthy"
################################################################################

covariate_data <- read.csv(SV_COVARIATE, sep = "\t")
covariate_data$sample.ID <- sub("^long", "", covariate_data$ID)

csv_files <- list.files(CSV_DECONVOLDED_PATH, pattern = "\\.csv$", full.names = TRUE)

# Read case legend and prepare group IDs
case_legend <- read.delim(CASE_LEGEND)
case_legend <- case_legend[case_legend$DiseaseCode %in% c("control", "FTD-C9"), ]

# 2. Renombrar valores dentro de DiseaseCode
case_legend$DiseaseCode[case_legend$DiseaseCode == "control"] <- "Healthy"
case_legend$DiseaseCode[case_legend$DiseaseCode == "FTD-C9"] <- "C9orf72"

case_legend$group_ID <- factor(case_legend$DiseaseCode)
case_legend$group.ID <- NULL

# Ensure levels are consistent and drop unused levels
conditions <- levels(as.factor(case_legend$group_ID))
conditions <- "C9orf72"

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
      comparison_data$group_ID <- factor(comparison_data$group_ID)
      
      comparison_data$SV1 <- covariate_data$SV1[match(comparison_data$RimodID, covariate_data$ID)]

      comparison_data$group_ID <- factor(comparison_data$group_ID, levels = c(REFERENCE_GROUP, condition))
      comparison_data <- comparison_data[order(comparison_data$group_ID), ]
      common_ids <- intersect(comparison_data$RimodID, rownames(data))
      data_subset <- data[rownames(data) %in% common_ids, ]
      data_ordered <- data_subset[match(comparison_data$RimodID, rownames(data_subset)), ]
      # CLEAN NA
      data_ordered[!is.finite(as.matrix(data_ordered))] <- 0
      data_ordered <- data_ordered[complete.cases(data_ordered), ]
      comparison_data <- comparison_data[match(rownames(data_ordered), comparison_data$RimodID), ]
      
      group <- factor(comparison_data$group_ID[comparison_data$group_ID %in% c(REFERENCE_GROUP, condition)])
      y <- DGEList(counts=t(data_ordered),group=group)
      y <- y[, colSums(y$counts) > 0]
      # Keep only samples in comparison_data that are in y
      comparison_data <- comparison_data[match(colnames(y), comparison_data$RimodID), ]
      
      # Average Log CPM histogram
      AveLogCPM <- aveLogCPM(y)
      y <- normLibSizes(y)
      
      # MD plot
      plot_filename <- file.path(condition_folder, "MD_plot_h.png")
      png(file = plot_filename)
      MD_plot <- plotMD(y, column=1) +
        abline(h=0, col="red", lty=2, lwd=2)
      dev.off()
      
      # DESIGN MATRIX
      group <- factor(comparison_data$group_ID, levels = c(REFERENCE_GROUP, condition))
      comparison_data$group_ID <- relevel(comparison_data$group_ID, ref = REFERENCE_GROUP)
      design <- model.matrix(~ 0 + comparison_data$group_ID + comparison_data$SV1, data = comparison_data)
      colnames(design) <- c(REFERENCE_GROUP, condition, "SV1")
      nrow(design) == ncol(y)
      
      # DISPERSION ESTIMATION
      y <- estimateDisp(y,design)
      fit <- glmFit(y,design)
      
      # Diferential expression analysis
      B.LvsP <- makeContrasts(paste0(condition, " - ", REFERENCE_GROUP), levels=design)
      res <- glmLRT(fit, contrast=B.LvsP)
      is.de <- decideTests(res, adjust.method = "fdr", p.value = 0.05, lfc= 0)
      results_table <- res$table
      results_table$adj_pval <- p.adjust(results_table$PValue, method = "fdr", n = length(results_table$PValue))
      results_table <- cbind(results_table, data.frame(is.de))
      # SAVING CSV RESULTS
      write.csv(res$table, file = file.path(condition_folder, "res_h.csv"))
      write.csv(results_table, file = file.path(condition_folder, "results_adj_h.csv"))
      ## Visualizations
      # "Volcano like" plot
      plot_filename <- file.path(condition_folder, "MD_res_plot_h.png")
      png(file = plot_filename)
      plotMD(res, status=is.de)
      dev.off()
      
      # Heatmap clustering 
      logCPM <- cpm(y, prior.count=2, log=TRUE)
      rownames(logCPM) <- y$genes$Symbol
      colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
      tr <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
      o <- order(tr$table$PValue)
      logCPM <- logCPM[o[1:30],]
      
      plot_filename <- file.path(condition_folder, "Heatmap_h.png")
      png(file = plot_filename)
      coolmap(logCPM, margins=c(7,7), lhei=c(1,6), lwid=c(1,3))
      dev.off()
      
      df[condition, tools::file_path_sans_ext(basename(file))] <- nrow(results_table[results_table$adj_pval<0.05,])
      
    }, error = function(e){
      df[condition, tools::file_path_sans_ext(basename(file))] <- 0
    })
  }
}

# Save summary DataFrame
write.csv(df, file = file.path(OUTPUT_DIRECTORY, "summary_h.csv"))


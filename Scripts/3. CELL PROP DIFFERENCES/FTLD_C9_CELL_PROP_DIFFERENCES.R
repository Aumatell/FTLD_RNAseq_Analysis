library("ggplot2")
library("openxlsx")
library("dplyr")
library("coin")  

#################################T.TESTS########################################
# FILE PATHS
frequency_file <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv"
metadata_file <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx"
output_path <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/FTLD/c9/"
################################################################################

if (!file.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
if (!file.exists(file.path(output_path,"Significative"))) {
  dir.create(file.path(output_path,"Significative"), recursive = TRUE)
}

data <- read.xlsx(metadata_file)

freq_cell_types <- read.csv(frequency_file)
freq_cell_types$X <- sapply(freq_cell_types$X, function(x) gsub("long", "", x))
freq_cell_types$X <- sapply(freq_cell_types$X, function(x) gsub("X", "", x))

merged_data <- merge(data, freq_cell_types, by.x = "sample.ID", by.y = "X")
merged_data[merged_data == 0] <- 1e-10

################################################################################
# Generate needed objects

m_conditions <- unique(merged_data$group.ID)
m_conditions_nh <- m_conditions[m_conditions != "Healthy"]
files_to_iterate <- colnames(merged_data)[4:length(colnames(merged_data))]
m_cond <- m_conditions_nh[1]
counter = 0

significatives <- data.frame(matrix(nrow = length(files_to_iterate), ncol = 3))
colnames(significatives)<- c("cell", "group+sv", "group")


for (column in colnames(freq_cell_types)[2:length(colnames(freq_cell_types))]){
  counter <- counter + 1
  filename <- paste0("Kruskal_", column, "_", m_cond, ".txt")
  
  subset_data <- merged_data[merged_data$group.ID %in% c("Healthy", m_cond), 
                             c("group.ID", column)]
  
  # Ensure group.ID is treated as a factor
  subset_data$group.ID <- factor(subset_data$group.ID, levels = c("Healthy", m_cond))
  
  # Calculate group means
  mean_healthy <- mean(subset_data[[column]][subset_data$group.ID == "Healthy"], na.rm = TRUE)
  mean_condition <- mean(subset_data[[column]][subset_data$group.ID == m_cond], na.rm = TRUE)
  
  # Calculate fold change
  fold_change <- ifelse(mean_healthy == 0, NA, mean_condition / mean_healthy)  # Avoid division by zero
  
  # Log-transform the fold change (optional)
  log2_fold_change <- ifelse(is.na(fold_change), NA, log2(fold_change))
  
  # Check if both groups (Healthy and m_cond) have at least 2 observations
  save = TRUE
  tryCatch({
    if (save == TRUE){sink(paste0(output_path, filename))}
    
    print(filename)
    
    # Model A: Kruskal-Wallis Test (without covariates)
    kruskal_A <- kruskal.test(subset_data[[column]], subset_data$group.ID)
    print("Kruskal-Wallis Test (General model)")
    print(kruskal_A)
    if (save == TRUE){sink()}
    
  }, error = function(e) {
    message("Error in processing column ", column, ": ", e$message)
    next
  })
  
  # Store the p-values, fold change, and means in the significatives data frame
  significatives[counter, 1] <- column
  significatives[counter, 3] <- kruskal_A$p.value  # Group (from Kruskal-Wallis test)
  significatives[counter, 4] <- fold_change  # Fold change
  significatives[counter, 5] <- log2_fold_change  # Log2 fold change
  significatives[counter, 6] <- mean_healthy  # Mean of the Healthy group
  significatives[counter, 7] <- mean_condition  # Mean of the specific condition group
}

# Update column names of the significatives data frame
colnames(significatives) <- c("cell", "group+sv", "group", "fold_change", 
                              "log2_fold_change", "mean_healthy", "mean_condition")

# Save the significatives data frame as a TSV file
tsv_output_path <- file.path(output_path, "significatives_nonparametric_with_means_cs.tsv")
write.table(significatives, 
            file = tsv_output_path, 
            sep = "\t",          
            row.names = FALSE,   
            col.names = TRUE,    
            quote = FALSE)

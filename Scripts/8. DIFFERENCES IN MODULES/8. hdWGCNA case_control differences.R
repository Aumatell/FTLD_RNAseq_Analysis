library(dplyr)
library(hdWGCNA)
library(xlsx)

for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", full.names = FALSE, recursive = FALSE)) {
  try({
    
    # Load Seurat object
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", CS, "/", CS, "_seurat.rds"))
    
    # Load metadata
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", 
                                row.names = 1, sheetIndex = 1)
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    # Harmonize sample names
    A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean")
    A@meta.data <- meta_joined
    
    # Get Module Eigengenes
    MEs <- hdWGCNA::GetMEs(A)
    
    # Ensure sample order matches between MEs and metadata
    sample_order <- rownames(MEs)
    # Ensure consistent factor order: Healthy (control), TDP (case)
    group_vector <- factor(
      A@meta.data[match(sample_order, A@meta.data$Sample), "group.ID"],
      levels = c("Healthy", "TDP")
    )
    
    # Apply Wilcoxon test and calculate fold change
    test_results <- lapply(1:ncol(MEs), function(i) {
      vec <- MEs[, i]
      group1 <- vec[group_vector == "Healthy"]
      group2 <- vec[group_vector == "TDP"]
      fc <- mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)  # TDP - Healthy
      test <- wilcox.test(group1, group2)
      
      list(
        module = colnames(MEs)[i],
        fold_change = fc,
        p.value = test$p.value
      )
    })
    
    cor_df <- do.call(rbind, lapply(test_results, as.data.frame))
    cor_df_sorted <- cor_df[order(cor_df$p.value), ]
    
    # Create output directory
    outdir <- paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", CS, "/MODULES/")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # Write gene lists for significant modules
    module_colors <- A@misc$pseudobulk$wgcna_degrees
    sig_modules <- cor_df_sorted$module[!is.na(cor_df_sorted$module)]
    
    for (mod in sig_modules) {
      genes <- module_colors$gene_name[module_colors$module == mod]
      write.csv(genes, file = paste0(outdir,mod, ".csv"), row.names = FALSE)
    }
    
    # Write summary and hub genes
    Hubs <- GetHubGenes(A, n_hubs = 10, mods = NULL, wgcna_name = NULL)
    write.csv(cor_df_sorted, file = paste0(outdir, "/cor_summary.csv"), row.names = FALSE)
    write.csv(as.data.frame(Hubs), file = paste0(outdir, "/hubs.csv"), row.names = FALSE)
    
    print(paste0("Finished: ", CS))
  })
}


library(Seurat)
library(dplyr)
library(hdWGCNA)

# Loop over directories
for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/NEW/CS_bo/", full.names = FALSE, recursive = FALSE)) {
  try({
    
    # Load Seurat object
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/NEW/CS_bo/", CS, "/", CS, "_seurat.rds"))
    
    # Load metadata
    metadata <- read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt", 
                           row.names = 1, sep = "\t")
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- gsub("-", "\\.",rownames(metadata_df))
    rownames(metadata_df) <- metadata_df$Sample_clean 
    # Filter out FTLD-TDP-C and recode group labels
    metadata_df <- metadata_df[metadata_df$GROUP != "FTLD-TDP-C", ]
    metadata_df$GROUP <- ifelse(metadata_df$GROUP == "Control", "Healthy", "TDP")
    
    samples_to_keep <- intersect(rownames(A@meta.data), metadata_df$Sample_clean)
    A <- subset(A,cells= samples_to_keep)
                
    # Harmonize and merge metadata
    A@meta.data <- metadata_df
    
    # Get Module Eigengenes
    MEs <- hdWGCNA::GetMEs(A)
    samples_in_A <- rownames(A@meta.data)
    MEs <- MEs[rownames(MEs) %in% samples_in_A, ]
    sample_order <- rownames(MEs)
    
    # Ensure correct group assignment order
    group_vector <- factor(
      A@meta.data[match(sample_order, rownames(A@meta.data)), "GROUP"],
      levels = c("Healthy", "TDP")
    )
    
    # Apply Wilcoxon test and compute fold change
    test_results <- lapply(1:ncol(MEs), function(i) {
      vec <- MEs[, i]
      group1 <- vec[group_vector == "Healthy"]
      group2 <- vec[group_vector == "TDP"]
      fc <- mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)  # TDP - Healthy
      test <- wilcox.test(group1, group2)
      
      list(
        module = colnames(MEs)[i],
        fold_change = fc,
        p.value = test$p.value
      )
    })
    
    # Compile and sort results
    cor_df <- do.call(rbind, lapply(test_results, as.data.frame))
    cor_df_sorted <- cor_df[order(cor_df$p.value), ]
    
    # Create output directory
    outdir <- paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/NEW/CS_bo/", CS, "/MODULES/")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # Write gene lists for significant modules
    module_colors <- A@misc$pseudobulk$wgcna_degrees
    sig_modules <- cor_df_sorted$module[!is.na(cor_df_sorted$module)]
    
    for (mod in sig_modules) {
      genes <- module_colors$gene_name[module_colors$module == mod]
      write.csv(genes, file = paste0(outdir, mod, ".csv"), row.names = FALSE)
    }
    
    # Write summary and hub genes
    Hubs <- GetHubGenes(A, n_hubs = 10, mods = NULL, wgcna_name = NULL)
    write.csv(cor_df_sorted, file = paste0(outdir, "/cor_summary.csv"), row.names = FALSE)
    write.csv(as.data.frame(Hubs), file = paste0(outdir, "/hubs.csv"), row.names = FALSE)
    
    print(paste0("Finished: ", CS))
  })
}

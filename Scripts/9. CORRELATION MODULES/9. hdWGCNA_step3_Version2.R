
############################## Correlation MEs VS Covariables ###############################
library(hdWGCNA)
library(dplyr)
library(Seurat)
library(xlsx)

CS <- list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", full.names = FALSE, recursive = FALSE)[1]

for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", 
                     full.names = FALSE, recursive = FALSE)) {
  if (file.exists(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", CS, "/", CS, "_seurat.rds"))) {
    
    tryCatch({
      
      message("Processing CS: ", CS)
#      sink(log_file, append = TRUE); cat(Sys.time(), " - Processing: ", CS, "\n"); sink()
      
      # Load Seurat object
      A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", CS, "/", CS, "_seurat.rds"))
      
      # Load metadata
      metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", 
                                  row.names = 1, sheetIndex = 1)
      covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/ArcSin/TDP/ArcSin_FTD_TDP_neuropath_SOM")
      
      # Merge metadata
      covariables$X <- gsub("long", "", covariables$X)
      A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
      metadata_df <- as.data.frame(metadata)
      metadata_df$Sample_clean <- rownames(metadata_df)
      
      meta_joined <- A@meta.data %>%
        mutate(Sample_clean = gsub("^X", "", Sample)) %>%
        left_join(metadata_df, by = "Sample_clean") %>%
        left_join(covariables, by = c("Sample_clean" = "X"))
      
      rownames(meta_joined) <- colnames(A)
      A@meta.data <- meta_joined

      # Align rownames
      if (!all(rownames(A@meta.data) == colnames(A))) {
        rownames(A@meta.data) <- colnames(A)
      }
      
      # Ensure WGCNA results are present
      if (!"pseudobulk" %in% names(A@misc)) stop("WGCNA results not found in @misc$pseudobulk for ", CS)
      
      # Remove outlier
      outlier <- "7BLACK"
      if (outlier %in% colnames(A)) {
        A <- subset(A, cells = setdiff(colnames(A), outlier))
        rownames(A@meta.data) <- colnames(A)
      }
      
      # Ensure numeric trait
      A@meta.data$TDP43b <- as.numeric(as.character(A@meta.data$TDP43b))
      if (all(is.na(A@meta.data$TDP43b))) stop("TDP43b column is all NA for ", CS)
      
      # Calculate module eigengenes inside Seurat object
      MEs <- hdWGCNA::GetMEs(A, wgcna_name = "pseudobulk")
      A@misc$pseudobulk$MEs <- MEs
      
      # Check MEs exist
      if (is.null(A@misc$pseudobulk$MEs)) stop("No MEs found in object for ", CS)
      
      # Run Module-Trait correlation
      cor_results <- hdWGCNA::ModuleTraitCorrelation(
        seurat_obj = A,
        traits = "TDP43b",
        wgcna_name = "pseudobulk"
      )
      
      # Extract correlation dataframe safely
      cor_df <- GetModuleTraitCorrelation(cor_results)
      
      # Sort correlations
      # Extract correlations for 'all_cells'
      cor_vec <- cor_df$cor$all_cells
      pval_vec <- cor_df$pval$all_cells
      fdr_vec <- cor_df$fdr$all_cells
      
      # Convert to data frame
      cor_df_clean <- data.frame(
        module = names(cor_vec),
        cor = as.numeric(cor_vec),
        p.value = as.numeric(pval_vec),
        fdr = as.numeric(fdr_vec)
      )
      
      # Sort by absolute correlation and then p-value
      library(dplyr)
      cor_df_sorted <- cor_df_clean %>%
        arrange(desc(abs(cor)), p.value)
      
      # Save results
      outdir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_asin_moduletraitcorrelation_tdp43b_tdp/"
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
      write.csv(cor_df_sorted, file = paste0(outdir, CS, ".csv"), row.names = FALSE)      
      # Cleanup
      rm(A, metadata, covariables, 
         cor_vec, pval_vec, fdr_vec,
         cor_df_sorted, cor_results, cor_df, cor_df_clean); gc()
      message("Finished processing CS: ", CS)
      
#      sink(log_file, append = TRUE); cat(Sys.time(), " - Finished: ", CS, "\n"); sink()
      
    }, error = function(e) {
#      sink(log_file, append = TRUE)
      cat(Sys.time(), "ERROR in ", CS, ": ", conditionMessage(e), "\n")
#      sink()
      message("Error in CS: ", CS, " → ", conditionMessage(e))
    })
  }
}



############################## Correlation MEs VS Covariables ###############################

# C9 

# Must change the covariable name for each case.

library(hdWGCNA)
library(dplyr)
library(Seurat)
library(xlsx)

CS <- list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/c9/CS_0.25_npcs2/", full.names = FALSE, recursive = FALSE)[1]

for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/c9/CS_0.25_npcs2/", 
                     full.names = FALSE, recursive = FALSE)) {
  if (file.exists(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/c9/CS_0.25_npcs2/", CS, "/", CS, "_seurat.rds"))) {
    
    tryCatch({
      
      message("Processing CS: ", CS)
      #      sink(log_file, append = TRUE); cat(Sys.time(), " - Processing: ", CS, "\n"); sink()
      
      # Load Seurat object
      A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/c9/CS_0.25_npcs2/", CS, "/", CS, "_seurat.rds"))
      
      # Load metadata
      metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", 
                                  row.names = 1, sheetIndex = 1)
      covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/ArcSin/C9/ArcSin_FTD_C9_neuropath_SOM.csv")
      
      # Merge metadata
      covariables$X <- gsub("X", "", covariables$X)
      A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
      metadata_df <- as.data.frame(metadata)
      metadata_df$Sample_clean <- rownames(metadata_df)
      
      meta_joined <- A@meta.data %>%
        mutate(Sample_clean = gsub("^X", "", Sample)) %>%
        left_join(metadata_df, by = "Sample_clean") %>%
        left_join(covariables, by = c("Sample_clean" = "X"))
      
      rownames(meta_joined) <- colnames(A)
      A@meta.data <- meta_joined
      
      # Align rownames
      if (!all(rownames(A@meta.data) == colnames(A))) {
        rownames(A@meta.data) <- colnames(A)
      }
      
      # Ensure WGCNA results are present
      if (!"pseudobulk" %in% names(A@misc)) stop("WGCNA results not found in @misc$pseudobulk for ", CS)
      
      # Remove outlier
      outlier <- ""
      if (outlier %in% colnames(A)) {
        A <- subset(A, cells = setdiff(colnames(A), outlier))
        rownames(A@meta.data) <- colnames(A)
      }
      
      # Ensure numeric trait
      # A@meta.data$TDP43b <- as.numeric(as.character(A@meta.data$TDP43b))
      # if (all(is.na(A@meta.data$TDP43b))) stop("TDP43b column is all NA for ", CS)
      
      # Calculate module eigengenes inside Seurat object
      MEs <- hdWGCNA::GetMEs(A, wgcna_name = "pseudobulk")
      A@misc$pseudobulk$MEs <- MEs
      
      # Check MEs exist
      if (is.null(A@misc$pseudobulk$MEs)) stop("No MEs found in object for ", CS)
      
      # Run Module-Trait correlation
      cor_results <- hdWGCNA::ModuleTraitCorrelation(
        seurat_obj = A,
        traits = "pTDP43",
        wgcna_name = "pseudobulk"
      )
      
      # Extract correlation dataframe safely
      cor_df <- GetModuleTraitCorrelation(cor_results)
      
      # Sort correlations
      # Extract correlations for 'all_cells'
      cor_vec <- cor_df$cor$all_cells
      pval_vec <- cor_df$pval$all_cells
      fdr_vec <- cor_df$fdr$all_cells
      
      # Convert to data frame
      cor_df_clean <- data.frame(
        module = names(cor_vec),
        cor = as.numeric(cor_vec),
        p.value = as.numeric(pval_vec),
        fdr = as.numeric(fdr_vec)
      )
      
      # Sort by absolute correlation and then p-value
      library(dplyr)
      cor_df_sorted <- cor_df_clean %>%
        arrange(desc(abs(cor)), p.value)
      
      # Save results
      outdir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/C9/CS_Moduletraitcorrelations_Asin_pTDP43/"
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
      write.csv(cor_df_sorted, file = paste0(outdir, CS, ".csv"), row.names = FALSE)      
      # Cleanup
      rm(A, metadata, covariables, 
         cor_vec, pval_vec, fdr_vec,
         cor_df_sorted, cor_results, cor_df, cor_df_clean); gc()
      message("Finished processing CS: ", CS)
      
      #      sink(log_file, append = TRUE); cat(Sys.time(), " - Finished: ", CS, "\n"); sink()
      
    }, error = function(e) {
      #      sink(log_file, append = TRUE)
      cat(Sys.time(), "ERROR in ", CS, ": ", conditionMessage(e), "\n")
      #      sink()
      message("Error in CS: ", CS, " → ", conditionMessage(e))
    })
  }
}

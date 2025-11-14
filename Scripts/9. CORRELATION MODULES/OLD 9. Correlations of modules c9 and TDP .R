library(dplyr)
library(xlsx)
library(Seurat)
library(hdWGCNA)
# TDP
all_correlations <- list()

# Recorrem els directoris
for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", 
                     full.names = FALSE, recursive = FALSE)) {
  try({
    message("Processant: ", CS)
    
    # Carreguem objecte Seurat
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", 
                        CS, "/", CS, "_seurat.rds"))
    
    # Metadata bulk i covariables
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", 
                                sheetIndex = 1)
    rownames(metadata) <- metadata[,1]
    metadata <- metadata[,-1]
    
    covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
    covariables$X <- gsub("long", "", covariables$X)
    
    # Selecció de mostres TDP
    selection <- rownames(metadata[metadata$group.ID == "TDP", , drop = FALSE])
    
    # Afegim columna Sample_clean
    A$Sample_clean <- gsub("X", "", A@meta.data$Sample)
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean") %>%
      left_join(covariables, by = c("Sample_clean" = "X"))
    
    # Afegim l’expressió dels gens d’interès
    expr <- FetchData(A, vars = c("STMN2", "TDP43b"))
    meta_joined <- cbind(meta_joined, expr)
    A@meta.data <- meta_joined
    
    # Eigengenes
    MEs <- hdWGCNA::GetMEs(A)
    MEs <- MEs[match(selection, gsub("X", "",rownames(MEs))), , drop = FALSE]
    
    # Correlacions
    CS_cor_list <- list()
    for (marker in c("STMN2", "TDP43b")) {
      if (!marker %in% colnames(meta_joined)) {
        warning("Marker ", marker, " no trobat a meta_joined per ", CS)
        next
      }
      y <- meta_joined[match(selection, gsub("X","",rownames(meta_joined))), marker]
      if (all(is.na(y))) {
        warning("Marker ", marker, " tot NA per ", CS)
        next
      }
      
      cor_results <- apply(MEs, 2, function(module) {
        cor.test(module, y)
      })
      
      cor_df <- data.frame(
        module = names(cor_results),
        cor = sapply(cor_results, function(x) x$estimate),
        p.value = sapply(cor_results, function(x) x$p.value),
        cell = CS,
        neuropath_marker = marker
      )
      CS_cor_list[[marker]] <- cor_df
    }
    
    # Combina i guarda
    if (length(CS_cor_list) > 0) {
      cor_df_combined <- bind_rows(CS_cor_list)
      
      outdir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_bo_ambhubs2/"
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      write.csv(cor_df_combined,
                file = paste0(outdir, CS, ".csv"),
                row.names = FALSE)
      
      all_correlations[[CS]] <- cor_df_combined
    }
  })
}

# Guardem totes les correlacions juntes
if (length(all_correlations) > 0) {
  all_cor_df <- bind_rows(all_correlations)
  write.csv(all_cor_df, "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/TDP_CS_combined_correlations.csv", row.names = FALSE)
}

# C9


all_correlations <- list()

for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/c9/CS_0.25_npcs2/", full.names = FALSE, recursive = FALSE)) {
  try({
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/c9/CS_0.25_npcs2/",CS,"/",CS,"_seurat.rds"))
    
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", row.names = 1, sheetIndex = 1)
    covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
    
    covariables$X <- gsub("X", "", covariables$X)
    A$Sample_clean <- gsub("X", "", A@meta.data$Sample)
    
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean") %>%
      left_join(covariables, by = c("Sample_clean" = "X"))
    
    A@meta.data <- meta_joined
    MEs <- hdWGCNA::GetMEs(A)
    
    # Store all marker correlations for this CS
    CS_cor_list <- list()
    
    for (marker in c("STMN2", "ACSL3", "lncRNA", "polyGP", "polyGA", "polyGR", "pTDP43", "fociSENSE", "fociANTI")) {
      cor_results <- apply(MEs, 2, function(module) {
        cor.test(module, as.numeric(meta_joined[, marker]))
      })
      
      cor_df <- data.frame(
        module = names(cor_results),
        cor = sapply(cor_results, function(x) x$estimate),
        p.value = sapply(cor_results, function(x) x$p.value),
        cell = CS,
        neuropath_marker = marker
      )
      
      CS_cor_list[[marker]] <- cor_df
    }
    
    # Combine all marker results for this CS
    cor_df_combined <- bind_rows(CS_cor_list)
    
    # Save per CS
    dir.create(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/C9/CS_bo/"), recursive = TRUE)
    write.csv(cor_df_combined,
              file = paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/C9/CS_bo/",CS,".csv"),
              row.names = FALSE)
    
    # Append to master list
    all_correlations[[CS]] <- cor_df_combined
  })
}

# Optional: Save all correlations across all CSs in one file
all_cor_df <- bind_rows(all_correlations)
write.csv(all_cor_df, "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/C9/C9_CS_combined_correlations.csv", row.names = FALSE)

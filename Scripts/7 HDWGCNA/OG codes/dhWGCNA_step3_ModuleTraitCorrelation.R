library(dplyr)
library(tidyverse)
library(hdWGCNA)
library(readxl)

# Define input/output
root_dir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/"
out_dir  <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_Final(test)/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (CS in list.dirs(root_dir, full.names = FALSE, recursive = FALSE)) {
  try({
    message("Processing: ", CS)
    seurat_path <- file.path(root_dir, CS, paste0(CS, "_seurat.rds"))
    if (!file.exists(seurat_path)) {
      message("No Seurat object found for ", CS)
      next
    }
    
    # Load Seurat object
    A <- readRDS(seurat_path)
    
    # Metadata
    metadata <- read_excel(
      "/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx",
      sheet = 1
    )
    metadata_df <- as.data.frame(metadata)
    rownames(metadata_df) <- metadata_df[,1]  # first col should be sample IDs
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    # Subset to only TDP samples
    metadata_df <- metadata_df[metadata_df$group.ID == "TDP", ]
    
    # Covariables
    covariables <- read.csv(
      "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv"
    )
    covariables$X <- gsub("long", "", covariables$X)
    
    # Harmonize and join
    A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean") %>%
      left_join(covariables, by = c("Sample_clean" = "X"))
    
    # Remove any rows that are still NA in key columns (like Sample_clean or TDP43b)
    meta_joined <- meta_joined[!is.na(meta_joined$Sample_clean) & !is.na(meta_joined$TDP43b), ]
    
    
    # Now restrict Seurat object to only TDP samples
    keep_cells <- meta_joined$Sample[meta_joined$group.ID == "TDP"]
    rownames(meta_joined) <- meta_joined$Sample
    A <- subset(A, cells = keep_cells)
    A@meta.data <- meta_joined[rownames(A@meta.data), , drop = FALSE]
    A@meta.data <- meta_joined[match(colnames(A), rownames(meta_joined)), , drop = FALSE]

    
    # Traits: select numeric columns to correlate (e.g. TDP43b)
    traits <- c("TDP43b", "STMN2")  
    A$STMN2 <- as.numeric(A$STMN2)
    A$TDP43b <- as.numeric(A$TDP43b)
    # Run ModuleTraitCorrelation only on TDP samples
    cor_df <- ModuleTraitCorrelation(
      seurat_obj = A,
      traits = traits,# "TDP43b",
      wgcna_name = "pseudobulk",   # adjust if you used another name
      cor_fun = "spearman",
      use = "pairwise.complete.obs"
    )
    
    cor_df_sorted <- cor_df %>% arrange(p.val)
    
    # Save results
    write.csv(cor_df_sorted,
              file = file.path(out_dir, paste0(CS, "_module_trait_cor_TDP.csv")),
              row.names = FALSE)
    
    # Optional: heatmap
    p_heatmap <- ModuleTraitHeatmap(
      seurat_obj = A,
      cor_df = cor_df,
      show_pval = TRUE
    )
    ggsave(
      filename = file.path(out_dir, paste0(CS, "_heatmap_TDP.png")),
      plot = p_heatmap,
      width = 8, height = 6, dpi = 300
    )
    
    message("Finished: ", CS)
  })
}

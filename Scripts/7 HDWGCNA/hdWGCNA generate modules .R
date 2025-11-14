# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# paralel processing
library(fs)
library(future.apply)

# Configurar sessio
plan(multisession, workers = parallel::detectCores() - 1)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

# Directori d'entrada i sortida
input_dir <- "/media/jaumatell/datos/URI/PROJECTE_SEURAT_BP/25_09/FC/Rimod_C9/CELL STATE/"
output_dir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/Rimod/CS/"

if (!dir_exists(output_dir)) dir_create(output_dir)
csv_files <- dir_ls(input_dir, glob = "*.csv")

WGCNA_pseudobulk <- function(file) {
  tryCatch({
    
    # Obtenir el nom base de l'arxiu sense extensió
    sample_name <- path_ext_remove(path_file(file))
    sample_output_dir <- file.path(output_dir, sample_name)
    
    if (!dir_exists(sample_output_dir)) dir_create(sample_output_dir)
    
    message("Processing sample: ", sample_name)
    
    # Carregar dades en format CSV
    sc_data <- read.csv(file, row.names = 1)
    
    # Convertir a objecte Seurat
    seurat_obj <- CreateSeuratObject(counts = t(sc_data))
    seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)
    
    seurat_obj@meta.data$Sample <- rownames(sc_data)
    seurat_obj@meta.data$cell_type <- sample_name
    
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      gene_select = "fraction",
      fraction = 0.05,
      wgcna_name = "pseudobulk"
    )
    
    message("Selected genes: ", length(GetWGCNAGenes(seurat_obj)))
    
    # log2CPM normalization
    cpm <- t(apply(sc_data, 1, function(x) {
      y <- x / sum(x) * 1e6
      log2(y + 1)
    }))
    
    seurat_obj <- SetDatExpr(seurat_obj, mat = cpm)
    seurat_obj@assays$RNA$data <- t(cpm)
    
    # Soft thresholding
    seurat_obj <- TestSoftPowers(seurat_obj)
    
    # Save power plot
    png(file.path(sample_output_dir, "soft_power.png"), width = 1000, height = 800)
    PlotSoftPowers(seurat_obj)
    dev.off()
    
    # Construct network
    seurat_obj <- ConstructNetwork(
      seurat_obj,
      tom_name = "pseudobulk",
      overwrite_tom = TRUE,
      mergeCutHeight = 0.25
    )
    
    # Dendrogram
    png(file.path(sample_output_dir, "dendrogram.png"), width = 1000, height = 800)
    PlotDendrogram(seurat_obj, main = "pseudobulk dendrogram")
    dev.off()
    
    # Eigengenes and connectivity
    seurat_obj <- ModuleEigengenes(seurat_obj, npcs = 2)
    seurat_obj <- ModuleConnectivity(seurat_obj)
    
    # DotPlot of MEs
    MEs <- GetMEs(seurat_obj)
    mods <- setdiff(colnames(MEs), "grey")
    meta <- seurat_obj@meta.data
    seurat_obj@meta.data <- cbind(meta, MEs)
    
    p_dot <- DotPlot(seurat_obj, features = mods, group.by = "cell_type") +
      RotatedAxis() +
      scale_color_gradient(high = "red", low = "grey95") +
      xlab("") + ylab("")
    
    png(file.path(sample_output_dir, "dotplot_MEs.png"), width = 1000, height = 800)
    print(p_dot)
    dev.off()
    
    # Reset metadata
    seurat_obj@meta.data <- meta
    
    # UMAP
    seurat_obj <- RunModuleUMAP(
      seurat_obj,
      n_hubs = 5,
      n_neighbors = 10,
      min_dist = 0.4,
      spread = 3,
      supervised = TRUE,
      target_weight = 0.3
    )
    
    umap_df <- GetModuleUMAP(seurat_obj)
    centroid_df <- umap_df %>%
      dplyr::group_by(module) %>%
      dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
    
    p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
      geom_point(color = umap_df$color, size = umap_df$kME * 2) +
      geom_label(data = centroid_df, label = as.character(centroid_df$module),
                 fontface = "bold", size = 2) +
      umap_theme() +
      theme(panel.background = element_rect(fill = "black"))
    
    png(file.path(sample_output_dir, "umap_modules.png"), width = 1000, height = 800)
    print(p_umap)
    dev.off()
    
    # Save the final Seurat object
    saveRDS(seurat_obj, file = file.path(sample_output_dir, paste0(sample_name, "_seurat.rds")))
    
    message("Finished processing: ", sample_name)
    
  }, error = function(e) {
    message("⚠️ Error in sample: ", file)
    message("Details: ", e$message)
  })
}




#future_lapply(csv_files, WGCNA_pseudobulk, future.seed = TRUE)
for (file in csv_files){
  WGCNA_pseudobulk(file)
}


library(dplyr)
library(hdWGCNA)
library(xlsx)

# Loop over directories
for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.25_bo/", full.names = FALSE, recursive = FALSE)) {
  try({
    
    # Load Seurat object
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.25_bo/", CS, "/", CS, "_seurat.rds"))
    
    # Load and prepare metadata
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", 
                                row.names = 1, sheetIndex = 1)
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    # Harmonize and merge metadata
    A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean")
    
    A@meta.data <- meta_joined
    
    # Get Module Eigengenes
    MEs <- hdWGCNA::GetMEs(A)
    
    # Correlation of MEs with group.ID
    cor_results <- apply(MEs, 2, function(module) {
      cor.test(module, as.numeric(as.factor(A@meta.data$group.ID)), 
               method = "spearman", exact = TRUE)
    })
    
    cor_df <- data.frame(
      module = names(cor_results),
      cor = sapply(cor_results, function(x) x$estimate),
      p.value = sapply(cor_results, function(x) x$p.value)
    )
    
    cor_df_sorted <- cor_df[order(cor_df$p.value), ]
    
    # Create output directory if not exists
    outdir <- paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.25_bo/", CS, "/MODULES/")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # Write gene lists for significant modules
    module_colors <- A@misc$pseudobulk$wgcna_degrees
    
    sig_modules <- cor_df_sorted$module[cor_df_sorted$p.value < 0.05]
    
    for (mod in sig_modules) {
      genes <- module_colors$gene_name[module_colors$module == mod]
      write.csv(genes, file = paste0(outdir, mod, ".csv"), row.names = FALSE)
    }
    write.csv(cor_df_sorted, file = paste0(outdir, "/cor_summary.csv"), row.names = FALSE)
    print(paste0("Finished: ", CS))
    
  }) # end try
}


library(dplyr)
library(tidyverse)

for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/", full.names = FALSE, recursive = FALSE)) {
  try({
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_bo_ambhubs/",CS,"/",CS,"_seurat.rds"))
    
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", row.names = 1, sheetIndex = 1)
    
    covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
    
    covariables$X<-gsub("long", "", covariables$X)
    A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
    
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean") %>% 
      left_join(covariables, by = c("Sample_clean" = "X"))
    
    A@meta.data <- meta_joined
    MEs <-hdWGCNA::GetMEs(A)
    
    # Correlación con Sample type to select top modules
    cor_results <- apply(MEs, 2, function(module) {
      cor.test(module, as.numeric(meta_joined$TDP43b))
    })
    cor_df <- data.frame(
      module = names(cor_results),
      cor = sapply(cor_results, function(x) x$estimate),
      p.value = sapply(cor_results, function(x) x$p.value)
    )
    
    cor_df_sorted <- cor_df[order(abs(cor_df$p.value), decreasing = F), ]
    cor_df_sorted
    
    dir.create(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_bo_ambhubs_2/"), recursive = T)
    write.csv(cor_df_sorted,
              file = paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_bo_ambhubs_2/",CS,".csv")
    ) 
    
  })
}



for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.2_npcs2/", full.names = FALSE, recursive = FALSE)) {
  try({
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.2_npcs2/",CS,"/",CS,"_seurat.rds"))
    
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", row.names = 1, sheetIndex = 1)
    
    covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
    
    covariables$X<-gsub("long", "", covariables$X)
    A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
    
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean") %>% 
      left_join(covariables, by = c("Sample_clean" = "X"))
    
    A@meta.data <- meta_joined
    MEs <-hdWGCNA::GetMEs(A)
    
    # Correlación con Sample type to select top modules
    cor_results <- apply(MEs, 2, function(module) {
      cor.test(module, as.numeric(meta_joined$TDP43b))
    })
    cor_df <- data.frame(
      module = names(cor_results),
      cor = sapply(cor_results, function(x) x$estimate),
      p.value = sapply(cor_results, function(x) x$p.value)
    )
    
    cor_df_sorted <- cor_df[order(abs(cor_df$p.value), decreasing = F), ]
    cor_df_sorted
    
    dir.create(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_0.2_npcs2/"), recursive = T)
    write.csv(cor_df_sorted,
              file = paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_0.2_npcs2/",CS,".csv")
    ) 
    
  })
}


for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.25_nou/", full.names = FALSE, recursive = FALSE)) {
  try({
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.25_nou/",CS,"/",CS,"_seurat.rds"))
    
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", row.names = 1, sheetIndex = 1)
    
    covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
    
    covariables$X<-gsub("long", "", covariables$X)
    A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
    
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean") %>% 
      left_join(covariables, by = c("Sample_clean" = "X"))
    
    A@meta.data <- meta_joined
    MEs <-hdWGCNA::GetMEs(A)
    
    # Correlación con Sample type to select top modules
    cor_results <- apply(MEs, 2, function(module) {
      cor.test(module, as.numeric(meta_joined$TDP43b))
    })
    cor_df <- data.frame(
      module = names(cor_results),
      cor = sapply(cor_results, function(x) x$estimate),
      p.value = sapply(cor_results, function(x) x$p.value)
    )
    
    cor_df_sorted <- cor_df[order(abs(cor_df$p.value), decreasing = F), ]
    cor_df_sorted
    
    dir.create(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_0.25_nou/"), recursive = T)
    write.csv(cor_df_sorted,
              file = paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_0.25_nou/",CS,".csv")
    ) 
    
  })
}



for (CS in list.dirs("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.5_nou/", full.names = FALSE, recursive = FALSE)) {
  try({
    A <- readRDS(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS_0.5_nou/",CS,"/",CS,"_seurat.rds"))
    
    metadata <- xlsx::read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", row.names = 1, sheetIndex = 1)
    
    covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
    
    covariables$X<-gsub("long", "", covariables$X)
    A$Sample_clean <- gsub("^X", "", A@meta.data$Sample)
    
    metadata_df <- as.data.frame(metadata)
    metadata_df$Sample_clean <- rownames(metadata_df)
    
    meta_joined <- A@meta.data %>%
      mutate(Sample_clean = gsub("^X", "", Sample)) %>%
      left_join(metadata_df, by = "Sample_clean") %>% 
      left_join(covariables, by = c("Sample_clean" = "X"))
    
    A@meta.data <- meta_joined
    MEs <-hdWGCNA::GetMEs(A)
    
    # Correlación con Sample type to select top modules
    cor_results <- apply(MEs, 2, function(module) {
      cor.test(module, as.numeric(meta_joined$TDP43b))
    })
    cor_df <- data.frame(
      module = names(cor_results),
      cor = sapply(cor_results, function(x) x$estimate),
      p.value = sapply(cor_results, function(x) x$p.value)
    )
    
    cor_df_sorted <- cor_df[order(abs(cor_df$p.value), decreasing = F), ]
    cor_df_sorted
    
    dir.create(paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_0.5_nou/"), recursive = T)
    write.csv(cor_df_sorted,
              file = paste0("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS_0.5_nou/",CS,".csv")
    ) 
    
  })
}
library(car)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(unikn)
library(cluster)
library(tidyverse)
library(readxl)

safe_cor <- function(x, y){
  tryCatch({
    if(all(is.na(x)) | all(is.na(y))) return(list(estimate = NA, p.value = NA))
    res <- cor.test(x, y, method = "spearman")
    list(estimate = res$estimate, p.value = res$p.value)
  }, error = function(e) list(estimate = NA, p.value = NA))
}


# C9
Expression_per_cell_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/CELL STATE ORIGINAL/"

OG_Covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
OG_Covariables_log    <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log/C9/Log_FTD_C9_neuropath_SOM.csv")
OG_Covariables_log2   <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log2/C9/Log2_FTD_C9_neuropath_SOM.csv")
OG_Covariables_logit  <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Logit/C9/Logit_FTD_C9_neuropath_SOM.csv")
OG_Covariables_arcsin <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/ArcSin/C9/ArcSin_FTD_C9_neuropath_SOM.csv")

# Excitatory neurons
cells <- c("CUX2_RASGRF2", "CUX2_RORB", "SCN4B_NEFH", "RORB_FOXO1", 
           "RORB_POU3F2", "RORB_ADGRL4", "RORB_LRRK1", "PCP4_NXPH2", 
           "VAT1L_EYA4", "VAT1L_THSD4", "THEMIS_NR4A2", "THEMIS_TMEM233", 
           "TLE4_CCBE1", "TLE4_MEGF11", "TLE4_SEMA3D")

genes <- c("STMN2", "NPTX2", "C9orf72", "CXCL10", "SPP1", "CX3CR1") 
covariables <- c("STMN2", "pTDP43")

results <- data.frame(Cell = character(),
                      Gene = character(),
                      Covariable = character(),
                      Gene_transformation = character(),
                      Cov_transformation = character(),
                      Spearman = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

for (cell in cells){
  data <- read.csv(paste0(Expression_per_cell_directory, cell, ".csv"), row.names = "X")
  
  for (gene in genes){
    common_ids <- intersect(rownames(data), OG_Covariables_log$X)
    subset_gene_subset <- data[common_ids, gene]
    subset_Covariables  <- OG_Covariables[match(common_ids, OG_Covariables$X), ]
    subset_Covariables_log  <- OG_Covariables_log[match(common_ids, OG_Covariables_log$X),]
    subset_Covariables_log2  <- OG_Covariables_log2[match(common_ids, OG_Covariables_log2$X), ]
    subset_Covariables_logit  <- OG_Covariables_logit[match(common_ids, OG_Covariables_logit$X), ]
    subset_Covariables_arcsin  <- OG_Covariables_arcsin[match(common_ids, OG_Covariables_arcsin$X), ]
    
    
    # gene transformations
    subset_log_data    <- log1p(gene_subset)
    subset_log2_data   <- log(gene_subset, base = 2)
    # protect logit and arcsin
    subset_logit_data <- tryCatch({
      logit(gene_subset, percents = FALSE, adjust = TRUE)
    }, error = function(e) rep(NA, length(gene_subset)))
    
    subset_norm_logit_data <- tryCatch({
      logit(gene_subset/max(gene_subset), percents = FALSE, adjust = TRUE)
    }, error = function(e) rep(NA, length(gene_subset)))
    
    subset_arcsin_data <- tryCatch({
      asin(sqrt(gene_subset))
    }, error = function(e) rep(NA, length(gene_subset)))
    
    subset_norm_arcsin_data <- tryCatch({
      asin(sqrt(gene_subset/max(gene_subset)))
    }, error = function(e) rep(NA, length(gene_subset)))
    
        
    for (covariable in covariables){
      
      if (covariable == "AAA"){
        
        outlier <- "7BLACK"
        no_outlier_ids <- common_ids[common_ids != "7BLACK"] 
        
        gene_subset <- subset_gene_subset[no_outlier_ids,]
        Covariables  <- subset_Covariables[match(no_outlier_ids, subset_Covariables$X), ]
        Covariables_log  <- subset_Covariables_log[match(no_outlier_ids, subset_Covariables_log$X),]
        Covariables_log2  <- subset_Covariables_log2[match(no_outlier_ids, subset_Covariables_log2$X), ]
        Covariables_logit  <- subset_Covariables_logit[match(no_outlier_ids, subset_Covariables_logit$X), ]
        Covariables_arcsin  <- subset_Covariables_arcsin[match(no_outlier_ids, subset_Covariables_arcsin$X), ]
        
        # gene transformations
        log_data    <- subset_log_data
        log2_data   <- subset_log2_data
        logit_data <- subset_logit_data
        norm_logit_data <- subset_norm_logit_data
        arcsin_data <- subset_arcsin_data
        norm_arcsin_data <- subset_norm_arcsin_data
        
      } else {
        gene_subset <- subset_gene_subset
        Covariables  <- subset_Covariables
        Covariables_log  <- subset_Covariables_log
        Covariables_log2  <- subset_Covariables_log2
        Covariables_logit  <- subset_Covariables_logit
        Covariables_arcsin  <- subset_Covariables_arcsin
        
        # gene transformations
        log_data    <- subset_log_data
        log2_data   <- subset_log2_data
        logit_data <- subset_logit_data
        norm_logit_data <- subset_norm_logit_data
        arcsin_data <- subset_arcsin_data
        norm_arcsin_data <- subset_norm_arcsin_data
      }
      
      result_corr <- safe_cor(gene_subset, Covariables[, covariable])
      
      result_log_corr <- safe_cor(log_data, Covariables[, covariable])
      result_log_log  <- safe_cor(log_data, Covariables_log[,covariable])
      result_gene_log <- safe_cor(gene_subset, Covariables_log[,covariable])
      
      result_log2_corr <- safe_cor(log2_data, Covariables[, covariable])
      result_log2_log2 <- safe_cor(log2_data, Covariables_log2[,covariable])
      result_gene_log2 <- safe_cor(gene_subset, Covariables_log2[,covariable])
      
      result_logit_corr <- safe_cor(logit_data, Covariables[, covariable])
      result_logit_logit <- safe_cor(logit_data, Covariables_logit[,covariable])
      result_gene_logit <- safe_cor(gene_subset, Covariables_logit[, covariable])
      
      result_norm_logit_corr <- safe_cor(norm_logit_data, Covariables[, covariable])
      result_norm_logit_logit <- safe_cor(norm_logit_data,  Covariables_logit[,covariable])
      
      result_arcsin_corr <- safe_cor(arcsin_data, Covariables[, covariable])
      result_arcsin_arcsin <- safe_cor(arcsin_data, Covariables_arcsin[,covariable])
      result_gene_arcsin <- safe_cor(gene_subset, Covariables_arcsin[,covariable])
      
      result_norm_arcsin_corr <- safe_cor(norm_arcsin_data, Covariables[, covariable])
      result_norm_arcsin_arcsin <- safe_cor(norm_arcsin_data, Covariables_arcsin[,covariable])
            
      results <- rbind(results,
                       #Cor - Cor
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "", 
                                  Spearman = result_corr$estimate, 
                                  p_value = result_corr$p.value),
                       #Log - Cor
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log",
                                  Cov_transformation = "", 
                                  Spearman = result_log_corr$estimate, 
                                  p_value = result_log_corr$p.value),
                       #Log - Log
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log",
                                  Cov_transformation = "Log", 
                                  Spearman = result_log_log$estimate, 
                                  p_value = result_log_log$p.value),
                       #Gene - Log
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Log", 
                                  Spearman = result_gene_log$estimate, 
                                  p_value = result_gene_log$p.value),
                       # Log2 - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log2",
                                  Cov_transformation = "", 
                                  Spearman = result_log2_corr$estimate, 
                                  p_value = result_log2_corr$p.value),
                       # Log2 - Log2
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log2",
                                  Cov_transformation = "Log2", 
                                  Spearman = result_log2_log2$estimate, 
                                  p_value = result_log2_log2$p.value),
                       # Gene - Log2
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Log2", 
                                  Spearman = result_gene_log2$estimate, 
                                  p_value = result_gene_log2$p.value),
                       # Logit - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Logit",
                                  Cov_transformation = "",  
                                  Spearman = result_logit_corr$estimate, 
                                  p_value = result_logit_corr$p.value),
                       # Logit - Logit
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Logit",
                                  Cov_transformation = "Logit",  
                                  Spearman = result_logit_logit$estimate, 
                                  p_value = result_logit_logit$p.value),
                       # Gene - Logit
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Logit",  
                                  Spearman = result_gene_logit$estimate, 
                                  p_value = result_gene_logit$p.value),
                       # Norm+Logit - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Logit",
                                  Cov_transformation = "",  
                                  Spearman = result_norm_logit_corr$estimate, 
                                  p_value = result_norm_logit_corr$p.value),
                       # Norm+Logit - Logit
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Logit",
                                  Cov_transformation = "Logit",  
                                  Spearman = result_norm_logit_logit$estimate, 
                                  p_value = result_norm_logit_logit$p.value),
                       # Arcsin - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Arcsin",
                                  Cov_transformation = "",  
                                  Spearman = result_arcsin_corr$estimate, 
                                  p_value = result_arcsin_corr$p.value),
                       # Arcsin - Arcsin
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Arcsin",
                                  Cov_transformation = "Arcsin",  
                                  Spearman = result_arcsin_arcsin$estimate, 
                                  p_value = result_arcsin_arcsin$p.value),
                       # Gene - Arcsin
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Arcsin",  
                                  Spearman = result_gene_arcsin$estimate, 
                                  p_value = result_gene_arcsin$p.value),
                       # Norm+Arcsin - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Arcsin",
                                  Cov_transformation = "",  
                                  Spearman = result_norm_arcsin_corr$estimate, 
                                  p_value = result_norm_arcsin_corr$p.value),
                       # Norm+Arcsin - Arcsin
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Arcsin",
                                  Cov_transformation = "Arcsin",  
                                  Spearman = result_norm_arcsin_arcsin$estimate, 
                                  p_value = result_norm_arcsin_arcsin$p.value))
      
      result_log_corr <-NULL
      result_log_log  <- NULL
      result_gene_log <- NULL
      result_log2_corr <- NULL
      result_log2_log2 <- NULL
      result_gene_log2 <- NULL
      result_logit_corr <- NULL
      result_logit_logit <- NULL
      result_gene_logit <- NULL
      result_norm_logit_corr <- NULL
      result_norm_logit_logit <- NULL
      result_arcsin_corr <- NULL
      result_arcsin_arcsin <- NULL
      result_gene_arcsin <- NULL
      result_norm_arcsin_corr <- NULL
      result_norm_arcsin_arcsin <- NULL
    }
  }
}

results


Covariables <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
Covariables_log    <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log/C9/Log_FTD_C9_neuropath_SOM.csv")
Covariables_log2   <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log2/C9/Log2_FTD_C9_neuropath_SOM.csv")
Covariables_logit  <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Logit/C9/Logit_FTD_C9_neuropath_SOM.csv")
Covariables_arcsin <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/ArcSin/C9/ArcSin_FTD_C9_neuropath_SOM.csv")

# Non-Neuronal
cells <- c("GFAP-neg", "GFAP-pos", "Micro", "Oligo", "OPC")
genes <- c("STMN2", "NPTX2", "C9orf72", "CXCL10", "SPP1", "CX3CR1")  
covariables <- c("STMN2", "pTDP43")


results <- data.frame(Cell = character(),
                      Gene = character(),
                      Covariable = character(),
                      Gene_transformation = character(),
                      Cov_transformation = character(),
                      Spearman = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

for (cell in cells){
  data <- read.csv(paste0(Expression_per_cell_directory, cell, ".csv"), row.names = "X")
  
  for (gene in genes){
    common_ids <- intersect(rownames(data), Covariables_log$X)
    gene_subset <- data[common_ids, gene]
    Covariables  <- Covariables[match(common_ids, Covariables$X), ]
    Covariables_log  <- Covariables_log[match(common_ids, Covariables_log$X),]
    Covariables_log2  <- Covariables_log2[match(common_ids, Covariables_log2$X), ]
    Covariables_logit  <- Covariables_logit[match(common_ids, Covariables_logit$X), ]
    Covariables_arcsin  <- Covariables_arcsin[match(common_ids, Covariables_arcsin$X), ]
    
    # Gene transformations
    
    # Log
    log_data    <- log1p(gene_subset)
    # Log2
    log2_data   <- log(gene_subset, base = 2)
    # Logit
    logit_data <- tryCatch({
      logit(gene_subset, percents = FALSE, adjust = TRUE)
    }, error = function(e) rep(NA, length(gene_subset)))
    # Norm+Logit
    norm_logit_data <- tryCatch({
      logit(gene_subset/max(gene_subset), percents = FALSE, adjust = TRUE)
    }, error = function(e) rep(NA, length(gene_subset)))
    # Arcsin
    arcsin_data <- tryCatch({
      asin(sqrt(gene_subset))
    }, error = function(e) rep(NA, length(gene_subset)))
    # Norm+Arcsin
    norm_arcsin_data <- tryCatch({
      asin(sqrt(gene_subset/max(gene_subset)))
    }, error = function(e) rep(NA, length(gene_subset)))
    
    for (covariable in covariables){
          common_ids <- intersect(rownames(data), OG_Covariables_log$X)
    subset_gene_subset <- data[common_ids, gene]
    subset_Covariables  <- OG_Covariables[match(common_ids, OG_Covariables$X), ]
    subset_Covariables_log  <- OG_Covariables_log[match(common_ids, OG_Covariables_log$X),]
    subset_Covariables_log2  <- OG_Covariables_log2[match(common_ids, OG_Covariables_log2$X), ]
    subset_Covariables_logit  <- OG_Covariables_logit[match(common_ids, OG_Covariables_logit$X), ]
    subset_Covariables_arcsin  <- OG_Covariables_arcsin[match(common_ids, OG_Covariables_arcsin$X), ]
    
    
    # gene transformations
    subset_log_data    <- log1p(gene_subset)
    subset_log2_data   <- log(gene_subset, base = 2)
    # protect logit and arcsin
    subset_logit_data <- tryCatch({
      logit(gene_subset, percents = FALSE, adjust = TRUE)
    }, error = function(e) rep(NA, length(gene_subset)))
    
    subset_norm_logit_data <- tryCatch({
      logit(gene_subset/max(gene_subset), percents = FALSE, adjust = TRUE)
    }, error = function(e) rep(NA, length(gene_subset)))
    
    subset_arcsin_data <- tryCatch({
      asin(sqrt(gene_subset))
    }, error = function(e) rep(NA, length(gene_subset)))
    
    subset_norm_arcsin_data <- tryCatch({
      asin(sqrt(gene_subset/max(gene_subset)))
    }, error = function(e) rep(NA, length(gene_subset)))
    
        
    for (covariable in covariables){
      
      if (covariable == "AAA"){
        
        outlier <- "7BLACK"
        no_outlier_ids <- common_ids[common_ids != "7BLACK"] 
        
        gene_subset <- subset_gene_subset[no_outlier_ids,]
        Covariables  <- subset_Covariables[match(no_outlier_ids, subset_Covariables$X), ]
        Covariables_log  <- subset_Covariables_log[match(no_outlier_ids, subset_Covariables_log$X),]
        Covariables_log2  <- subset_Covariables_log2[match(no_outlier_ids, subset_Covariables_log2$X), ]
        Covariables_logit  <- subset_Covariables_logit[match(no_outlier_ids, subset_Covariables_logit$X), ]
        Covariables_arcsin  <- subset_Covariables_arcsin[match(no_outlier_ids, subset_Covariables_arcsin$X), ]
        
        # gene transformations
        log_data    <- subset_log_data
        log2_data   <- subset_log2_data
        logit_data <- subset_logit_data
        norm_logit_data <- subset_norm_logit_data
        arcsin_data <- subset_arcsin_data
        norm_arcsin_data <- subset_norm_arcsin_data
        
      } else {
        gene_subset <- subset_gene_subset
        Covariables  <- subset_Covariables
        Covariables_log  <- subset_Covariables_log
        Covariables_log2  <- subset_Covariables_log2
        Covariables_logit  <- subset_Covariables_logit
        Covariables_arcsin  <- subset_Covariables_arcsin
        
        # gene transformations
        log_data    <- subset_log_data
        log2_data   <- subset_log2_data
        logit_data <- subset_logit_data
        norm_logit_data <- subset_norm_logit_data
        arcsin_data <- subset_arcsin_data
        norm_arcsin_data <- subset_norm_arcsin_data
      }
      
      result_corr <- safe_cor(gene_subset, Covariables[, covariable])
      
      result_log_corr <- safe_cor(log_data, Covariables[, covariable])
      result_log_log  <- safe_cor(log_data, Covariables_log[,covariable])
      result_gene_log <- safe_cor(gene_subset, Covariables_log[,covariable])
      
      result_log2_corr <- safe_cor(log2_data, Covariables[, covariable])
      result_log2_log2 <- safe_cor(log2_data, Covariables_log2[,covariable])
      result_gene_log2 <- safe_cor(gene_subset, Covariables_log2[,covariable])
      
      result_logit_corr <- safe_cor(logit_data, Covariables[, covariable])
      result_logit_logit <- safe_cor(logit_data, Covariables_logit[,covariable])
      result_gene_logit <- safe_cor(gene_subset, Covariables_logit[, covariable])
      
      result_norm_logit_corr <- safe_cor(norm_logit_data, Covariables[, covariable])
      result_norm_logit_logit <- safe_cor(norm_logit_data,  Covariables_logit[,covariable])
      
      result_arcsin_corr <- safe_cor(arcsin_data, Covariables[, covariable])
      result_arcsin_arcsin <- safe_cor(arcsin_data, Covariables_arcsin[,covariable])
      result_gene_arcsin <- safe_cor(gene_subset, Covariables_arcsin[,covariable])
      
      result_norm_arcsin_corr <- safe_cor(norm_arcsin_data, Covariables[, covariable])
      result_norm_arcsin_arcsin <- safe_cor(norm_arcsin_data, Covariables_arcsin[,covariable])
      
      results <- rbind(results,
                       #Cor - Cor
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "", 
                                  Spearman = result_corr$estimate, 
                                  p_value = result_corr$p.value),
                       #Log - Cor
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log",
                                  Cov_transformation = "", 
                                  Spearman = result_log_corr$estimate, 
                                  p_value = result_log_corr$p.value),
                       #Log - Log
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log",
                                  Cov_transformation = "Log", 
                                  Spearman = result_log_log$estimate, 
                                  p_value = result_log_log$p.value),
                       #Gene - Log
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Log", 
                                  Spearman = result_gene_log$estimate, 
                                  p_value = result_gene_log$p.value),
                       # Log2 - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log2",
                                  Cov_transformation = "", 
                                  Spearman = result_log2_corr$estimate, 
                                  p_value = result_log2_corr$p.value),
                       # Log2 - Log2
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Log2",
                                  Cov_transformation = "Log2", 
                                  Spearman = result_log2_log2$estimate, 
                                  p_value = result_log2_log2$p.value),
                       # Gene - Log2
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Log2", 
                                  Spearman = result_gene_log2$estimate, 
                                  p_value = result_gene_log2$p.value),
                       # Logit - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Logit",
                                  Cov_transformation = "",  
                                  Spearman = result_logit_corr$estimate, 
                                  p_value = result_logit_corr$p.value),
                       # Logit - Logit
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Logit",
                                  Cov_transformation = "Logit",  
                                  Spearman = result_logit_logit$estimate, 
                                  p_value = result_logit_logit$p.value),
                       # Gene - Logit
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Logit",  
                                  Spearman = result_gene_logit$estimate, 
                                  p_value = result_gene_logit$p.value),
                       # Norm+Logit - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Logit",
                                  Cov_transformation = "",  
                                  Spearman = result_norm_logit_corr$estimate, 
                                  p_value = result_norm_logit_corr$p.value),
                       # Norm+Logit - Logit
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Logit",
                                  Cov_transformation = "Logit",  
                                  Spearman = result_norm_logit_logit$estimate, 
                                  p_value = result_norm_logit_logit$p.value),
                       # Arcsin - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Arcsin",
                                  Cov_transformation = "",  
                                  Spearman = result_arcsin_corr$estimate, 
                                  p_value = result_arcsin_corr$p.value),
                       # Arcsin - Arcsin
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Arcsin",
                                  Cov_transformation = "Arcsin",  
                                  Spearman = result_arcsin_arcsin$estimate, 
                                  p_value = result_arcsin_arcsin$p.value),
                       # Gene - Arcsin
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "",
                                  Cov_transformation = "Arcsin",  
                                  Spearman = result_gene_arcsin$estimate, 
                                  p_value = result_gene_arcsin$p.value),
                       # Norm+Arcsin - Corr
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Arcsin",
                                  Cov_transformation = "",  
                                  Spearman = result_norm_arcsin_corr$estimate, 
                                  p_value = result_norm_arcsin_corr$p.value),
                       # Norm+Arcsin - Arcsin
                       data.frame(Cell = cell, Gene = gene, 
                                  Covariable = covariable,
                                  Gene_transformation = "Norm + Arcsin",
                                  Cov_transformation = "Arcsin",  
                                  Spearman = result_norm_arcsin_arcsin$estimate, 
                                  p_value = result_norm_arcsin_arcsin$p.value))
      
      result_log_corr <-NULL
      result_log_log  <- NULL
      result_gene_log <- NULL
      result_log2_corr <- NULL
      result_log2_log2 <- NULL
      result_gene_log2 <- NULL
      result_logit_corr <- NULL
      result_logit_logit <- NULL
      result_gene_logit <- NULL
      result_norm_logit_corr <- NULL
      result_norm_logit_logit <- NULL
      result_arcsin_corr <- NULL
      result_arcsin_arcsin <- NULL
      result_gene_arcsin <- NULL
      result_norm_arcsin_corr <- NULL
      result_norm_arcsin_arcsin <- NULL
      
    }
  }
  }
}
View(results)
results


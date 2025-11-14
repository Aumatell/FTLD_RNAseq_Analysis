
library("Seurat")
library("dplyr")
library("RColorBrewer")
library("unikn")
library("cluster")
library("tidyverse")
library("readxl")  # Asegurar que readxl está cargado

# Load Metadata
metadata <- read_xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx")
colnames(metadata)[1] <- "X"
samples <- metadata$X[metadata$group.ID == "C9orf72"]
metadata <- metadata[metadata$group.ID == "C9orf72",]

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
CRscores$X <- gsub(pattern = "X", replacement = "", x = CRscores$X)
CRscores <- CRscores[CRscores$X %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), CRscores$X)
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% filter(X %in% common_samples) %>% column_to_rownames("X")

if (nrow(proportions_data) > 2 && nrow(covariables) > 2) {
  
  results <- expand.grid(Cell_State = colnames(proportions_data), 
                         Covariate = colnames(covariables), 
                         stringsAsFactors = FALSE)
  
  compute_spearman <- function(cell_state, covariate) {
    test_result <- cor.test(proportions_data[, cell_state], covariables[, covariate], 
                            method = "spearman", exact = TRUE)
    return(c(test_result$estimate, test_result$p.value))
  }
  
  cor_pvals <- mapply(compute_spearman, results$Cell_State, results$Covariate)
  
  results$Spearman_Correlation <- cor_pvals[1, ]
  results$p_value <- cor_pvals[2, ]
}

write.csv(results, "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/C9_spearman_results_true.csv", row.names = FALSE)


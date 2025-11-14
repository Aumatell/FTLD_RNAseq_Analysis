library("Seurat")
library("dplyr")
library("RColorBrewer")
library("unikn")
library("cluster")
library("tidyverse")
library("readxl")  
library("car")

# Function to calculate correlation
compute_spearman <- function(cell_state, covariate, props, covars) {
  test_result <- cor.test(props[, cell_state], covars[, covariate],
                          method = "spearman", exact = TRUE)
  return(data.frame(Cell_State = cell_state,
                    Covariate = covariate,
                    Spearman_Correlation = test_result$estimate,
                    p_value = test_result$p.value))
}


# Load Metadata
metadata <- read_xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx")
colnames(metadata)[1] <- "X"
samples <- metadata$X[metadata$group.ID == "TDP"]
metadata <- metadata[metadata$group.ID == "TDP",]

# LOG
outdir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/TRANSFORMED_CELL_PROPORTION_CORRELATIONS/"

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log/TDP/Log_FTD_TDP_neuropath_SOM.csv")
CRscores$X <- gsub(pattern = "long", replacement = "", x = CRscores$X)
CRscores <- CRscores[CRscores$X %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/Log/FTLD/TDP/Log_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), CRscores$X)
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(X %in% common_samples) %>% column_to_rownames("X")

# Separate STMN2 to remove outlier
cov_stmn2 <- covariables[rownames(covariables) != "7BLACK", "STMN2", drop = FALSE]
proportions_data_STMN2 <- proportions_data[rownames(proportions_data) != "7BLACK", ]


cov_tdp43b <- covariables[, "TDP43b", drop = FALSE]

# Results dataframe 
results <- expand.grid(Cell_State = colnames(proportions_data), 
                       Covariate = colnames(covariables), 
                       stringsAsFactors = FALSE)

# --- Correlations with TDP43b ---
res_tdp43b <- do.call(rbind, lapply(colnames(proportions_data), function(cs) {
  compute_spearman(cs, "TDP43b", proportions_data, cov_tdp43b)
}))


# --- Correlations with STMN2 (without outlier) ---
res_stmn2 <- do.call(rbind, lapply(colnames(proportions_data_STMN2), function(cs) {
  compute_spearman(cs, "STMN2", proportions_data_STMN2, cov_stmn2)
}))

# Combine results
results <- rbind(res_tdp43b, res_stmn2)

write.csv(results, paste0(outdir, "Log_TDP_spearman_results.csv"), row.names = FALSE)


# LOG2

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log2/TDP/Log2_FTD_TDP_neuropath_SOM")
CRscores$X <- gsub(pattern = "long", replacement = "", x = CRscores$X)
CRscores <- CRscores[CRscores$X %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/Log2/FTLD/TDP/Log2_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), CRscores$X)
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(X %in% common_samples) %>% column_to_rownames("X")

# ULL AMB AQUEST CANVI DE VALOR
y <- covariables[, "STMN2"]
y[is.infinite(y)] <- min(y[is.finite(y)]) - 1e-2
covariables[, "STMN2"] <- y

# Separate STMN2 to remove outlier
cov_stmn2 <- covariables[rownames(covariables) != "7BLACK", "STMN2", drop = FALSE]
proportions_data_STMN2 <- proportions_data[rownames(proportions_data) != "7BLACK", ]


cov_tdp43b <- covariables[, "TDP43b", drop = FALSE]

# Results dataframe 
results <- expand.grid(Cell_State = colnames(proportions_data), 
                       Covariate = colnames(covariables), 
                       stringsAsFactors = FALSE)

# --- Correlations with TDP43b ---
res_tdp43b <- do.call(rbind, lapply(colnames(proportions_data), function(cs) {
  compute_spearman(cs, "TDP43b", proportions_data, cov_tdp43b)
}))


# --- Correlations with STMN2 (without outlier) ---
res_stmn2 <- do.call(rbind, lapply(colnames(proportions_data_STMN2), function(cs) {
  compute_spearman(cs, "STMN2", proportions_data_STMN2, cov_stmn2)
}))

# Combine results
results <- rbind(res_tdp43b, res_stmn2)

write.csv(results, paste0(outdir, "Log2_TDP_spearman_results.csv"), row.names = FALSE)


# LOGIT

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Logit/TDP/Logit_FTD_TDP_neuropath_SOM.csv")
CRscores$X <- gsub(pattern = "long", replacement = "", x = CRscores$X)
CRscores <- CRscores[CRscores$X %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/Logit/FTLD/TDP/Logit_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), CRscores$X)
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(X %in% common_samples) %>% column_to_rownames("X")

# ULL AMB AQUEST CANVI DE VALOR
y <- covariables[, "STMN2"]
y[is.infinite(y)] <- min(y[is.finite(y)]) - 1e-2
covariables[, "STMN2"] <- y

# Separate STMN2 to remove outlier
cov_stmn2 <- covariables[rownames(covariables) != "7BLACK", "STMN2", drop = FALSE]
proportions_data_STMN2 <- proportions_data[rownames(proportions_data) != "7BLACK", ]


cov_tdp43b <- covariables[, "TDP43b", drop = FALSE]

# Results dataframe 
results <- expand.grid(Cell_State = colnames(proportions_data), 
                       Covariate = colnames(covariables), 
                       stringsAsFactors = FALSE)

# --- Correlations with TDP43b ---
res_tdp43b <- do.call(rbind, lapply(colnames(proportions_data), function(cs) {
  compute_spearman(cs, "TDP43b", proportions_data, cov_tdp43b)
}))


# --- Correlations with STMN2 (without outlier) ---
res_stmn2 <- do.call(rbind, lapply(colnames(proportions_data_STMN2), function(cs) {
  compute_spearman(cs, "STMN2", proportions_data_STMN2, cov_stmn2)
}))

# Combine results
results <- rbind(res_tdp43b, res_stmn2)


write.csv(results, paste0(outdir, "Logit_TDP_spearman_results.csv"), row.names = FALSE)



# ARCSIN

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/ArcSin/TDP/ArcSin_FTD_TDP_neuropath_SOM")
CRscores$X <- gsub(pattern = "long", replacement = "", x = CRscores$X)
CRscores <- CRscores[CRscores$X %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/ArcSin/FTLD/TDP/Arcsin_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), CRscores$X)
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(X %in% common_samples) %>% column_to_rownames("X")

# ULL AMB AQUEST CANVI DE VALOR
y <- covariables[, "STMN2"]
y[is.infinite(y)] <- min(y[is.finite(y)]) - 1e-2
covariables[, "STMN2"] <- y

# Separate STMN2 to remove outlier
cov_stmn2 <- covariables[rownames(covariables) != "7BLACK", "STMN2", drop = FALSE]
proportions_data_STMN2 <- proportions_data[rownames(proportions_data) != "7BLACK", ]


cov_tdp43b <- covariables[, "TDP43b", drop = FALSE]

# Results dataframe 
results <- expand.grid(Cell_State = colnames(proportions_data), 
                       Covariate = colnames(covariables), 
                       stringsAsFactors = FALSE)

# --- Correlations with TDP43b ---
res_tdp43b <- do.call(rbind, lapply(colnames(proportions_data), function(cs) {
  compute_spearman(cs, "TDP43b", proportions_data, cov_tdp43b)
}))


# --- Correlations with STMN2 (without outlier) ---
res_stmn2 <- do.call(rbind, lapply(colnames(proportions_data_STMN2), function(cs) {
  compute_spearman(cs, "STMN2", proportions_data_STMN2, cov_stmn2)
}))

# Combine results
results <- rbind(res_tdp43b, res_stmn2)

write.csv(results, paste0(outdir, "ArcSin_TDP_spearman_results.csv"), row.names = FALSE)



# C9


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

# LOG
outdir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/TRANSFORMED_CELL_PROPORTION_CORRELATIONS/"

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log/C9/Log_FTD_C9_neuropath_SOM.csv", row.names = "X")
rownames(CRscores) <- gsub(pattern = "X", replacement = "", x = rownames(CRscores))
CRscores <- CRscores[rownames(CRscores) %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/Log/FTLD/C9/Log_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), rownames(CRscores))
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(rownames(CRscores) %in% common_samples) 


if (nrow(proportions_data) > 2 && nrow(covariables) > 2) {
  
  results <- expand.grid(Cell_State = colnames(proportions_data), 
                         Covariate = colnames(covariables), 
                         stringsAsFactors = FALSE)
  
  compute_spearman <- function(cell_state, covariate) {
    test_result <- cor.test(proportions_data[, cell_state], covariables[, covariate], 
                            method = "spearman", exact = TRUE) 
    return(c(test_result$estimate, test_result$p.value))
  }
  
  cor_pvals <- apply(results, 1, function(row) {
    compute_spearman(row["Cell_State"], row["Covariate"])
  })
  
  results$Spearman_Correlation <- cor_pvals[1, ]
  results$p_value <- cor_pvals[2, ]
}
write.csv(results, paste0(outdir, "Log_C9_spearman_results.csv"), row.names = FALSE)


# LOG2

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Log2/C9/Log2_FTD_C9_neuropath_SOM.csv", row.names = "X")
rownames(CRscores) <- gsub(pattern = "X", replacement = "", x = rownames(CRscores))
CRscores <- CRscores[rownames(CRscores) %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/Log2/FTLD/C9/Log2_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), rownames(CRscores))
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(rownames(CRscores) %in% common_samples) 
if (nrow(proportions_data) > 2 && nrow(covariables) > 2) {
  
  results <- expand.grid(Cell_State = colnames(proportions_data), 
                         Covariate = colnames(covariables), 
                         stringsAsFactors = FALSE)
  
  compute_spearman <- function(cell_state, covariate) {
    test_result <- cor.test(proportions_data[, cell_state], covariables[, covariate], 
                            method = "spearman", exact = TRUE) 
    return(c(test_result$estimate, test_result$p.value))
  }
  
  cor_pvals <- apply(results, 1, function(row) {
    compute_spearman(row["Cell_State"], row["Covariate"])
  })
  
  results$Spearman_Correlation <- cor_pvals[1, ]
  results$p_value <- cor_pvals[2, ]
}

write.csv(results, paste0(outdir, "Log2_C9_spearman_results.csv"), row.names = FALSE)


# LOGIT

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/Logit/C9/Logit_FTD_C9_neuropath_SOM.csv", row.names = "X")
rownames(CRscores) <- gsub(pattern = "X", replacement = "", x = rownames(CRscores))
CRscores <- CRscores[rownames(CRscores) %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/Logit/FTLD/C9/Logit_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), rownames(CRscores))
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(rownames(CRscores) %in% common_samples) 
if (nrow(proportions_data) > 2 && nrow(covariables) > 2) {
  
  results <- expand.grid(Cell_State = colnames(proportions_data), 
                         Covariate = colnames(covariables), 
                         stringsAsFactors = FALSE)
  
  compute_spearman <- function(cell_state, covariate) {
    test_result <- cor.test(proportions_data[, cell_state], covariables[, covariate], 
                            method = "spearman", exact = TRUE) 
    return(c(test_result$estimate, test_result$p.value))
  }
  
  cor_pvals <- apply(results, 1, function(row) {
    compute_spearman(row["Cell_State"], row["Covariate"])
  })
  
  results$Spearman_Correlation <- cor_pvals[1, ]
  results$p_value <- cor_pvals[2, ]
}
write.csv(results, paste0(outdir, "Logit_C9_spearman_results.csv"), row.names = FALSE)



# ARCSIN

CRscores <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/TRANSFORMED_COVARIABLES/ArcSin/C9/ArcSin_FTD_C9_neuropath_SOM.csv", row.names = "X")
rownames(CRscores) <- gsub(pattern = "X", replacement = "", x = rownames(CRscores))
CRscores <- CRscores[rownames(CRscores) %in% samples,]

proportions_directory <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/TRANSFORMATED_PROPORTIONS/ArcSin/FTLD/C9/Arcsin_Cell_state_proportions.csv"
proportions_data <- read.csv(proportions_directory, row.names = 1)
rownames(proportions_data) <- gsub("^X", "", rownames(proportions_data))

################################################################################
common_samples <- intersect(rownames(proportions_data), rownames(CRscores))
proportions_data <- proportions_data[common_samples, , drop = FALSE]
covariables <- CRscores %>% dplyr::filter(rownames(CRscores) %in% common_samples) 
if (nrow(proportions_data) > 2 && nrow(covariables) > 2) {
  
  results <- expand.grid(Cell_State = colnames(proportions_data), 
                         Covariate = colnames(covariables), 
                         stringsAsFactors = FALSE)
  
  compute_spearman <- function(cell_state, covariate) {
    test_result <- cor.test(proportions_data[, cell_state], covariables[, covariate], 
                            method = "spearman", exact = TRUE) 
    return(c(test_result$estimate, test_result$p.value))
  }
  
  cor_pvals <- apply(results, 1, function(row) {
    compute_spearman(row["Cell_State"], row["Covariate"])
  })
  
  results$Spearman_Correlation <- cor_pvals[1, ]
  results$p_value <- cor_pvals[2, ]
}

write.csv(results, paste0(outdir, "ArcSin_C9_spearman_results.csv"), row.names = FALSE)


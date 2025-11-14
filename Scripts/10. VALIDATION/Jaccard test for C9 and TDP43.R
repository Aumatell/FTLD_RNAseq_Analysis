# JACCARD 
library("openxlsx")
library("GeneOverlap")
library("VennDiagram")
library("grid")

# Define functions
# Jaccard similarity index
jaccard_similarity <- function(list1, list2) {
  intersection_size <- length(intersect(list1, list2))
  union_size <- length(union(list1, list2))
  return(intersection_size / union_size)
}

# Perform permutation test for p-value
permutation_test <- function(list1, list2, n_permutations = 10000) {
  observed_jaccard <- jaccard_similarity(list1, list2)
  
  permuted_jaccards <- replicate(n_permutations, {
    permuted_list1 <- sample(list1)
    permuted_list2 <- sample(list2)
    jaccard_similarity(permuted_list1, permuted_list2)
  })
  p_value <- mean(permuted_jaccards >= observed_jaccard)
  return(p_value)
}

hypergeometric_test <- function(list1, list2, total_population_size) {
  intersection_size <- length(intersect(list1, list2))
  total_genes <- total_population_size  # Total number of possible genes
  
  # Number of DEGs in list1 and list2
  K1 <- length(list1)
  K2 <- length(list2)
  
  # Hypergeometric test to calculate p-value
  p_value <- phyper(intersection_size - 1, K1, total_genes - K1, K2, lower.tail = FALSE)
  
  return(p_value)
}


# TDP43

# Define variables
CStates <- c(
  "Arterial", "DISC1_RELN", "Pericyte", "SMC", "TLE4_CCBE1",
  "Capillary", "GFAP-neg", "PVALB_CEMIP", "SST_ADAMTS19", "TLE4_MEGF11",
  "CDH4_CCK", "GFAP-pos", "PVALB_MYBPC1", "SST_BRINP3", "TLE4_SEMA3D",
  "CDH4_SCGN", "LAMP5_CA3", "PVALB_PTHLH", "SST_GALNT14", "VAT1L_EYA4",
  "CLMP_KCNMA1", "LAMP5_PMEPA1", "RORB_ADGRL4", "SST_NPY", "Venous",
  "CLMP_PDGFRA", "Micro", "RORB_FOXO1", "VIP_CLSTN2",
  "CUX2_RASGRF2", "Oligo", "RORB_LRRK1", "T_Cell", "VIP_HTR2C",
  "CUX2_RORB", "OPC", "RORB_POU3F2", "THEMIS_NR4A2", "VIP_LAMA3",
  "DISC1_CCK", "PCP4_NXPH2", "SCN4B_NEFH", "THEMIS_TMEM233"
)

FTLD_TDP <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/FTLD/TDP/CS_LRT/"
NOU <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/NEW/CS_SVcorrect/"
RESULTS_PATH <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/OVERLAP/TDP_NEW/CS_Correcte"

# Generate output dataframe
Results <- data.frame(matrix(ncol = length(CStates), nrow = 10))
colnames(Results) <- CStates
rownames(Results) <- c(
  "Upregulated FTLD_TDP", "Upregulated Nou", 
  "Downregulated FTLD_TDP", "Downregulated Nou",
  "CU FTLD – Nou", "CD FTLD – Nou", 
  "JU FTLD – Nou", "JD FTLD – Nou",
  "GU FTLD – Nou", "GD FTLD – Nou")

# Iterate to obtain the results for each case.
for (i in seq_along(CStates)) {
  CS <- CStates[i]
    # Read data
  if (file.exists(paste0(FTLD_TDP, CS ,"/", "TDP/results_adj_h.csv"))){
    ftld_tdp <- read.csv(paste0(FTLD_TDP, CS ,"/", "TDP/results_adj_h.csv"))
  } else {
    ftld_tdp <- c()
  }
  
  if (file.exists(paste0(NOU,CS,"/" ,"TDP/results_adj_h.csv"))){
    nou <- read.csv(paste0(NOU,CS,"/" ,"TDP/results_adj_h.csv"))
  } else {
    nou <- c()
  }
  
  # Extract DEG

  # For U_ftld
  if (is.null(ftld_tdp$X[ftld_tdp$X1.TDP..1.Healthy == 1])) {
    U_ftld <- character(0)  # Assign an empty vector if NULL
  } else {
    U_ftld <- ftld_tdp$X[ftld_tdp$X1.TDP..1.Healthy == 1]
  }
  
  # For U_nou
  if (is.null(nou$X[nou$X.1.Control.1.TDP == 1])) {
    U_nou <- character(0)  # Assign an empty vector if NULL
  } else {
    U_nou <- nou$X[nou$X.1.Control.1.TDP == 1]
  }
  
  # For D_ftld
  if (is.null(ftld_tdp$X[ftld_tdp$X1.TDP..1.Healthy == -1])) {
    D_ftld <-character(0)  # Assign an empty vector if NULL
  } else {
    D_ftld <- ftld_tdp$X[ftld_tdp$X1.TDP..1.Healthy == -1]
  }
  
  # For D_nou
  if (is.null(nou$X[nou$X.1.Control.1.TDP == -1])) {
    D_nou <-character(0) # Assign an empty vector if NULL
  } else {
    D_nou <- nou$X[nou$X.1.Control.1.TDP == -1]
  }
  
  # OVERLAP GENES
  C_U_ftld_nou <- intersect(U_ftld, U_nou)
  C_D_ftld_nou <- intersect(D_ftld, D_nou)

  # SIMILARITY INDEX
  Ujaccard13 <- jaccard_similarity(U_ftld, U_nou)
  Djaccard13 <- jaccard_similarity(D_ftld, D_nou)

  # HIPERGEOMETRIC TESTS
  total_genes_FN <- length(unique(c(ftld_tdp$X, nou$X)))

  # Gene Overlap
  GU_13 = testGeneOverlap(newGeneOverlap(U_ftld, U_nou, total_genes_FN))@pval
  GD_13 = testGeneOverlap(newGeneOverlap(D_ftld, D_nou, total_genes_FN))@pval

  # UPDATE RESULTS TABLE WITH THE CALCULATED VALUES.
  Results["Upregulated FTLD_TDP", CS] <- length(U_ftld)
  Results["Upregulated Nou", CS] <- length(U_nou)
  Results["Downregulated FTLD_TDP", CS] <- length(D_ftld)
  Results["Downregulated Nou", CS] <- length(D_nou)
  Results["CU FTLD – Nou", CS] <- length(C_U_ftld_nou)
  Results["CD FTLD – Nou", CS] <- length(C_D_ftld_nou)
  Results["JU FTLD – Nou", CS] <- Ujaccard13
  Results["JD FTLD – Nou", CS] <- Djaccard13
  Results["GU FTLD – Nou", CS] <- GU_13
  Results["GD FTLD – Nou", CS] <- GD_13
  
}
Results <- t(Results)
View(Results)
View(Results[rowSums(Results[, 1:4]) > 0, ])

write.csv(Results, file = paste0(RESULTS_PATH, "/TDP43_Validation_results_cs.csv"))



# C9

FTLD_TDP <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/FTLD/C9/CS_SV/"
NOU <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/RIMOD/CS_SV/"
RESULTS_PATH <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/OVERLAP/C9_NEW/CS"

# Generate output dataframe
Results <- data.frame(matrix(ncol = length(CStates), nrow = 10))
colnames(Results) <- CStates
rownames(Results) <- c(
  "Upregulated FTLD_TDP", "Upregulated Rimod", 
  "Downregulated FTLD_TDP", "Downregulated Rimod",
  "CU FTLD – Rimod", "CD FTLD – Rimod", 
  "JU FTLD – Rimod", "JD FTLD – Rimod",
  "GU FTLD – Rimod", "GD FTLD – Rimod")

# Iterate to obtain the results for each case.
for (i in seq_along(CStates)) {
  CS <- CStates[i]
  # Read data
  if (file.exists(paste0(FTLD_TDP, CS ,"/", "C9orf72/results_adj_h.csv"))){
    ftld_tdp <- read.csv(paste0(FTLD_TDP, CS ,"/", "C9orf72/results_adj_h.csv"))
  } else {
    ftld_tdp <- c()
  }
  
  if (file.exists(paste0(NOU,CS,"/" ,"C9orf72/results_adj_h.csv"))){
    nou <- read.csv(paste0(NOU,CS,"/" ,"C9orf72/results_adj_h.csv"))
  } else {
    nou <- c()
  }
  
  # Extract DEG
  
  # For U_ftld
  if (is.null(ftld_tdp$X[ftld_tdp$X.1.Healthy.1.C9orf72 == 1])) {
    U_ftld <- character(0)  # Assign an empty vector if NULL
  } else {
    U_ftld <- ftld_tdp$X[ftld_tdp$X.1.Healthy.1.C9orf72 == 1]
  }
  
  # For U_nou
  if (is.null(nou$X[nou$X.1.Healthy.1.C9orf72 == 1])) {
    U_nou <- character(0)  # Assign an empty vector if NULL
  } else {
    U_nou <- nou$X[nou$X.1.Healthy.1.C9orf72 == 1]
  }
  
  # For D_ftld
  if (is.null(ftld_tdp$X[ftld_tdp$X.1.Healthy.1.C9orf72 == -1])) {
    D_ftld <-character(0)  # Assign an empty vector if NULL
  } else {
    D_ftld <- ftld_tdp$X[ftld_tdp$X.1.Healthy.1.C9orf72 == -1]
  }
  
  # For D_nou
  if (is.null(nou$X[nou$X.1.Healthy.1.c9orf72 == -1])) {
    D_nou <-character(0) # Assign an empty vector if NULL
  } else {
    D_nou <- nou$X[nou$X.1.Healthy.1.C9orf72 == -1]
  }
  
  # OVERLAP GENES
  C_U_ftld_nou <- intersect(U_ftld, U_nou)
  C_D_ftld_nou <- intersect(D_ftld, D_nou)
  
  # SIMILARITY INDEX
  Ujaccard13 <- jaccard_similarity(U_ftld, U_nou)
  Djaccard13 <- jaccard_similarity(D_ftld, D_nou)
  
  # HIPERGEOMETRIC TESTS
  total_genes_FN <- length(unique(c(ftld_tdp$X, nou$X)))
  
  # Gene Overlap
  GU_13 = testGeneOverlap(newGeneOverlap(U_ftld, U_nou, total_genes_FN))@pval
  GD_13 = testGeneOverlap(newGeneOverlap(D_ftld, D_nou, total_genes_FN))@pval
  
  # UPDATE RESULTS TABLE WITH THE CALCULATED VALUES.
  Results["Upregulated FTLD_TDP", CS] <- length(U_ftld)
  Results["Upregulated Rimod", CS] <- length(U_nou)
  Results["Downregulated FTLD_TDP", CS] <- length(D_ftld)
  Results["Downregulated Rimod", CS] <- length(D_nou)
  Results["CU FTLD – Rimod", CS] <- length(C_U_ftld_nou)
  Results["CD FTLD – Rimod", CS] <- length(C_D_ftld_nou)
  Results["JU FTLD – Rimod", CS] <- Ujaccard13
  Results["JD FTLD – Rimod", CS] <- Djaccard13
  Results["GU FTLD – Rimod", CS] <- GU_13
  Results["GD FTLD – Rimod", CS] <- GD_13
  
}
Results <- t(Results)
View(Results)
View(Results[rowSums(Results[, 1:4]) > 0, ])
write.csv(Results, file = paste0(RESULTS_PATH, "/C9_Validation_results_cs.csv"))

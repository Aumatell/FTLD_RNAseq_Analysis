# GENES P-VAl

Gene_list <- c("TUBA4A","NEFL","UBQLN2","OPTN","CHCHD10","MATR3","TBK1",
               "C9orf72","VCP","CHMP2B","SQSTM1","TARDBP","HNRNPA1",
               "FUS","MAPT","DPP6","TMEM106B","GRN","ARPP21","UNC13A",
               "C19orf52","FARP2","TINAG","MZT1","TNIP1","RCL1",
               "PDS5B","C3AR1","SMG8","VIPR1","L3MBTL1","RBPJL",
               "ANO9","HNRNPL", "NPTX2")

directory_C9 <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/FTLD/C9/CS_LRT"
directory_TDP <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/FTLD/TDP/CS_LRT"

# Function to extract folder name (last folder before filename)
extract_foldername <- function(path) {
  parts <- strsplit(path, "/")[[1]]
  # Second to last element is folder
  folder <- parts[length(parts)-2]
  return(folder)
  }  
  
# C9
list_files <- list.files(directory_C9, pattern = "results_adj_h", full.names = TRUE, recursive = TRUE)

# Initialize results with gene column
results <- data.frame(gene = Gene_list, stringsAsFactors = FALSE)



# Loop over files
for (file in list_files) {
  celltypename <- extract_foldername(file)
  data <- read.csv(file, row.names = 1)
  
  # Match gene list to rownames of data
  subset <- data[match(Gene_list, rownames(data)), "PValue"]
  
  # Add to results
  results[[celltypename]] <- subset
}

results$Count_below_0.05 <- apply(results[,-1], 1, function(x) sum(x < 0.05, na.rm = TRUE))


write.csv(results, file = "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/C9_Pval_gens_URI.csv")

# TDP
list_files <- list.files(directory_TDP, pattern = "results_adj_h", full.names = TRUE, recursive = TRUE)

# Initialize results with gene column
results <- data.frame(gene = Gene_list, stringsAsFactors = FALSE)

# Function to extract folder name (last folder before filename)
extract_foldername <- function(path) {
  parts <- strsplit(path, "/")[[1]]
  # Second to last element is folder
  folder <- parts[length(parts)-2]
  return(folder)
}

# Loop over files
for (file in list_files) {
  celltypename <- extract_foldername(file)
  data <- read.csv(file, row.names = 1)
  
  # Match gene list to rownames of data
  subset <- data[match(Gene_list, rownames(data)), "PValue"]
  
  # Add to results
  results[[celltypename]] <- subset
}
results$Count_below_0.05 <- apply(results[,-1], 1, function(x) sum(x < 0.05, na.rm = TRUE))


write.csv(results, file = "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/TDP_Pval_gens_URI.csv")


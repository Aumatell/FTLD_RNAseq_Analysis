A <-read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CELL_PROP/FTLD/TDP/significatives_nonparametric_with_means_cs.tsv")



# Defineix el directori on tens els fitxers
directori <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS/TDP43"  

# Llista de fitxers (assumint que són .tsv o .csv, adapta segons calgui)
fitxers <- list.files(path = directori, pattern = "\\.txt$|\\.tsv$|\\.csv$", full.names = TRUE)

# Inicialitza una llista per guardar els resultats
resultats <- data.frame(
  fitxer = character(),
  pvalue_significatius = integer(),
  total_files = integer(),
  stringsAsFactors = FALSE
)

# Itera sobre cada fitxer
for (fitxer in fitxers) {
  # Llegeix el fitxer (ajusta el separador si cal)
  taula <- read.csv(fitxer, header = TRUE)
  
  # Compta les files amb p.value < 0.05
  significatius <- sum(taula$p.value < 0.05, na.rm = TRUE)
  
  # Compta el total de files
  total <- nrow(taula)
  
  # Afegeix a la taula de resultats
  resultats <- rbind(resultats, data.frame(
    fitxer = basename(fitxer),
    pvalue_significatius = significatius,
    total_files = total,
    stringsAsFactors = FALSE
  ))
}

# Mostra el resultat
print(resultats)
write.csv(resultats, file = "/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/EGG_VS_COV/FTLD/TDP/CS/resum_resultats_TDP43.csv", row.names = FALSE)



# Defineix el directori on hi ha les carpetes de cell states
directori_base <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/CS"  # Canvia això per la teva ruta

# Llista de subcarpetes (cell states)
cell_states <- list.dirs(path = directori_base, recursive = FALSE, full.names = TRUE)

# Inicialitza la taula de resultats
resultats <- data.frame(
  cell_state = character(),
  pvalue_significatius = integer(),
  total_moduls = integer(),
  stringsAsFactors = FALSE
)

# Itera per cada carpeta de cell state
for (carpeta in cell_states) {
  fitxer <- file.path(carpeta, "MODULES/cor_summary.csv")
  
  if (file.exists(fitxer)) {
    # Llegeix l’arxiu
    taula <- read.csv(fitxer, header = TRUE)
    
    # Compta mòduls significatius
    significatius <- sum(taula$p.value < 0.05, na.rm = TRUE)
    
    # Compta mòduls totals
    total <- nrow(taula)
    
    # Nom del cell state (nom de la carpeta)
    cell_state_nom <- basename(carpeta)
    
    # Afegeix-ho a la taula de resultats
    resultats <- rbind(resultats, data.frame(
      cell_state = cell_state_nom,
      pvalue_significatius = significatius,
      total_moduls = total,
      stringsAsFactors = FALSE
    ))
  }
}

# Guarda el resultat en un CSV
write.csv(resultats, file = "/media/jaumatell/datos/URI/BAYESPRISM_12_3/hdWGCNA/FTLD/TDP/resum_correlacio_moduls_Healthy_vs_ctrl_CS.csv", row.names = FALSE)


# Defineix el directori base amb les carpetes de tipus cel·lular
directori_base <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/NEW/CS"  # Substitueix amb la ruta correcta

# Llista de subcarpetes (tipus cel·lular)
carpetes <- list.dirs(path = directori_base, recursive = FALSE, full.names = TRUE)

# Inicialitza la taula de resultats
resultats <- data.frame(
  tipus_cellular = character(),
  n_positius = integer(),
  n_negatius = integer(),
  stringsAsFactors = FALSE
)

# Itera per cada carpeta
for (carpeta in carpetes) {
  fitxer <- file.path(carpeta, "FTLD/results_adj_h.csv")
  
  if (file.exists(fitxer)) {
    # Llegeix el fitxer
    taula <- read.csv(fitxer, header = TRUE)
    
    # Compta valors 1 i -1 en la columna "X1.FTLD..1.Control"
    positius <- sum(taula$X1.FTLD..1.Control == 1, na.rm = TRUE)
    negatius <- sum(taula$X1.FTLD..1.Control == -1, na.rm = TRUE)
    
    # Nom del tipus cel·lular (nom de la carpeta)
    nom_cell <- basename(carpeta)
    
    # Afegeix a la taula de resultats
    resultats <- rbind(resultats, data.frame(
      tipus_cellular = nom_cell,
      n_positius = positius,
      n_negatius = negatius,
      stringsAsFactors = FALSE
    ))
  }
}

# Guarda el resultat en un fitxer CSV
write.csv(resultats, file = "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/NEW/CS/resum_resultats_up_down.csv", row.names = FALSE)


###### VOLCANOPLOTS
library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)

results_dir <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/EDGER/FTLD/TDP/CS"
output_dir <- file.path(results_dir, "volcano_plots")
dir.create(output_dir, showWarnings = FALSE)

files <- list.files(results_dir, pattern = "results_adj_h\\.csv$", full.names = TRUE, recursive = TRUE)

for (file in files) {
  file_base <- basename(dirname(file))
  df <- tryCatch({
    read_csv(file, show_col_types = FALSE)
  }, error = function(e) {
    warning(paste("Error reading", file, ":", e$message))
    return(NULL)
  })
  if (is.null(df)) next
  colnames(df) <- c("Gene", "logFC", "logCPM", "LR", "PValue", "adj_pval", "TDP_vs_Healthy")
  required_cols <- c("logFC", "PValue", "adj_pval", "Gene")
  if (!all(required_cols %in% colnames(df))) {
    warning(paste("Skipping", file, "- missing required columns"))
    next
  }
  df$Gene <- as.character(df$Gene)
  p <- EnhancedVolcano(df,
                       lab = df$Gene,
                       x = 'logFC',
                       y = 'PValue',
                       title = paste("Volcano Plot -", file_base),
                       pCutoff = 0.05,
                       FCcutoff = 0,
                       pointSize = 2.5,
                       labSize = 3.5,
                       selectLab = df$Gene)
  ggsave(filename = file.path(output_dir, paste0(file_base, "_volcano.png")),
         plot = p, width = 8, height = 6, dpi = 300)
}

####################### VENN DIAGRAM OVERLAP 


# Install necessary packages if not already installed
packages <- c("eulerr", "patchwork", "purrr")
install_if_missing <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(install_if_missing)) install.packages(install_if_missing)

# Load libraries
library(eulerr)
library(patchwork)
library(purrr)

# Sample data frame with group-wise counts
predf <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/OVERLAP/TDP/CS/Validation_results_cs_FTLD_NEW.csv")
df <- data.frame(
  Group = predf$X,
  Healthy = predf$Upregulated.Nou,           # Total Healthy in each group
  FTLD_TDP = predf$Upregulated.FTLD_TDP,           # Total FTLD-TDP in each group
  Common = predf$CU.FTLD...Nou              # Overlap between Healthy and FTLD-TDP
)

# Generate Euler diagrams for each group
plots <- pmap(df, function(Group, Healthy, FTLD_TDP, Common) {
  fit <- euler(c(
    "Pottier" = Healthy - Common,
    "Sant Pau" = FTLD_TDP - Common,
    "Common genes" = Common
  ))
  
  plot(fit,
       main = Group,
       fills = list(fill = c("#1b9e77", "#d95f02"), alpha = 0.6),
       quantities = TRUE)
})

# Combine all plots using patchwork (like facet_wrap)
wrap_plots(plots) + plot_annotation(title = "Faceted Venn Diagrams by Group")


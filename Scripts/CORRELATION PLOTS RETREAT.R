

##### CELL PROP DIF PLOT
######################################################################
## 
# metadatas
SP_metadata <- readxl::read_xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx")
case_legend <- read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt", sep =  "\t", row.names = 1)
case_legend <- case_legend[case_legend$GROUP != "FTLD-TDP-C",]
case_legend$GROUP[case_legend$GROUP != "Control"] <- "TDP"

# RORB_LRRK1 (Sant pau)
sp1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/theta.state_original.csv", row.names = 1)

# RORB_LRRK1 (Validació)
n1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/theta.state_cellstate.csv", row.names = 1)

# NPTX2 (Sant pau)
sp2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)

# NPTX2 (Validació)
n2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)




# Match metadata and sample order
sp1_ids <- gsub("X", "", rownames(sp1))
SP_metadata_filtered <- SP_metadata[SP_metadata$sample.ID %in% sp1_ids, ]
SP_metadata_filtered <- SP_metadata_filtered[match(sp1_ids, SP_metadata_filtered$sample.ID), ]
Pheno <- SP_metadata_filtered$group.ID
Pheno[Pheno == "Healthy"] <- "Control"

# List of target cell states to visualize
target_cell_states <- c("GFAP.pos", "PVALB_PTHLH", "LAMP5_PMEPA1", "PVALB_CEMIP")
nice_names <- c("GFAP+", "PVALB PTHLH", "LAMP5 PMEPA1", "PVALB CEMIP")

# Initialize dataframe
cell_data <- data.frame()

# Loop through each cell type and collect data
for (i in seq_along(target_cell_states)) {
  cell_type <- target_cell_states[i]
  label <- nice_names[i]
  cell_prop <- sp1[rownames(sp1) %in% paste0("X", SP_metadata_filtered$sample.ID), cell_type]
  
  temp_df <- data.frame(
    Phenotype = Pheno,
    Cell_proportion = cell_prop,
    Cell_state = label
  )
  
  cell_data <- rbind(cell_data, temp_df)
}

# Ensure Phenotype is factor and ordered
cell_data$Phenotype <- factor(cell_data$Phenotype, levels = c("Control", "TDP"))

library(dplyr)
library(ggplot2)

# Step 1: Compute max y per facet (cell state)
max_y_per_state <- cell_data %>%
  group_by(Cell_state) %>%
  summarise(y = max(Cell_proportion, na.rm = TRUE) * 0.95)

# Step 2: Create annotation data frame
annotation_df <- data.frame(
  Cell_state = c("GFAP+", "PVALB PTHLH", "LAMP5 PMEPA1", "PVALB CEMIP"),
  label = c(
    "P = 0.003\nFC = 7.62",
    "P = 0.017\nFC = 0.72",
    "P = 0.023\nFC = 0.59",
    "P = 0.030\nFC = 0.62"
  ),
  x = 1.5  # halfway between Control and TDP
)

# Step 3: Merge y positions into annotation df
annotation_df <- left_join(annotation_df, max_y_per_state, by = "Cell_state")

# Step 4: Plot
ggplot(cell_data, aes(x = Phenotype, y = Cell_proportion, fill = Phenotype)) +
  geom_boxplot(size = 0.4, outlier.alpha = 0.3, alpha = 0.5) +
  facet_wrap(~Cell_state, scales = "free", nrow = 1) +
  labs(x = "", y = "Cell type proportion (%)") +
  scale_fill_manual(values = c("Control" = "deepskyblue4", "TDP" = "darkgoldenrod2")) +
  geom_text(data = annotation_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.5) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    panel.spacing = unit(5, "pt"),
    legend.position = "right"
  )


# C9
# metadatas
SP_metadata <- readxl::read_xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx")
case_legend <- read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt", sep =  "\t", row.names = 1)
case_legend <- case_legend[case_legend$GROUP != "FTLD-TDP-C",]
case_legend$GROUP[case_legend$GROUP != "Control"] <- "TDP"

# RORB_LRRK1 (Sant pau)
sp1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)

# RORB_LRRK1 (Validació)
n1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/theta.state_cellstate.csv", row.names = 1)

# NPTX2 (Sant pau)
sp2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)

# NPTX2 (Validació)
n2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)




# Match metadata and sample order
sp1_ids <- gsub("X", "", rownames(sp1))
SP_metadata_filtered <- SP_metadata[SP_metadata$sample.ID %in% sp1_ids, ]
SP_metadata_filtered <- SP_metadata_filtered[match(sp1_ids, SP_metadata_filtered$sample.ID), ]
Pheno <- SP_metadata_filtered$group.ID
Pheno[Pheno == "Healthy"] <- "Control"

# List of target cell states to visualize
target_cell_states <- c("VAT1L_EYA4", "GFAP.pos", "T_Cell", "PVALB_MYBPC1", "RORB_POU3F2")
nice_names <- c("VAT1L EYA4", "GFAP+", "T Cell", "PVALB MYBPC1", "RORB POU3F2")

# Initialize dataframe
cell_data <- data.frame()

# Loop through each cell type and collect data
for (i in seq_along(target_cell_states)) {
  cell_type <- target_cell_states[i]
  label <- nice_names[i]
  cell_prop <- sp1[rownames(sp1) %in% paste0("X", SP_metadata_filtered$sample.ID), cell_type]
  
  temp_df <- data.frame(
    Phenotype = Pheno,
    Cell_proportion = cell_prop,
    Cell_state = label
  )
  
  cell_data <- rbind(cell_data, temp_df)
}

# Ensure Phenotype is factor and ordered
cell_data$Phenotype <- factor(cell_data$Phenotype, levels = c("Control", "C9orf72"))

library(dplyr)
library(ggplot2)

# Step 1: Compute max y per facet (cell state)
max_y_per_state <- cell_data %>%
  group_by(Cell_state) %>%
  summarise(y = max(Cell_proportion, na.rm = TRUE) * 0.95)

# Step 2: Create annotation data frame
annotation_df <- data.frame(
  Cell_state = c("VAT1L EYA4", "GFAP+", "T Cell", "PVALB MYBPC1", "RORB POU3F2"),
  label = c(
    "P = 0.013\nFC = 0.067",
    "P = 0.026\nFC = 202.85",
    "P = 0.033\nFC = 1.24",
    "P = 0.009\nFC = 1.49",
    "P = 0.042\nFC = 0.87"
  ),
  x = 1.5  # halfway between Control and C9
)

# Step 3: Merge y positions into annotation df
annotation_df <- left_join(annotation_df, max_y_per_state, by = "Cell_state")

# Step 4: Plot
ggplot(cell_data, aes(x = Phenotype, y = Cell_proportion, fill = Phenotype)) +
  geom_boxplot(size = 0.4, outlier.alpha = 0.3, alpha = 0.5) +
  facet_wrap(~Cell_state, scales = "free", nrow = 1) +
  labs(x = "", y = "Cell type proportion (%)") +
  scale_fill_manual(values = c("Control" = "deepskyblue4", "C9orf72" = "darkgoldenrod2")) +
  geom_text(data = annotation_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.5) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    panel.spacing = unit(5, "pt"),
    legend.position = "right"
  )


############################################################################
#
# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/theta.state_original.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  GFAP_pos = A1$GFAP.pos,
  LAMP5_PMEPA1 = A1$LAMP5_PMEPA1,
  RORB_LRRK1 = A1$RORB_LRRK1,
  pTDP43 = C$TDP43b
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(GFAP_pos, LAMP5_PMEPA1, RORB_LRRK1),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(pTDP43, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(pTDP43, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.00007  # adjust as needed
corr_stats$y <- c(0.23, 0.18, 0.13)  # one y per cell type for visual spacing

# Define inverted colors
inverted_colors <- c(
  "GFAP_pos" = "#00BFC4",     # cyan
  "LAMP5_PMEPA1" = "green", 
  "RORB_LRRK1" = "#F8766D"    # red
)

# Plot with annotations
ggplot(SP_long, aes(x = pTDP43, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs pTDP-43 density",
    x = "pTDP-43 density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0, )



# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Assume `A` and `C` already read in as before
# A: cell proportions
# C: covariates with STMN2

C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

# Build the data
SP_stmn2 <- data.frame(
  TLE4_CCBE1 = A1$TLE4_CCBE1,
  CUX2_RORB = A1$CUX2_RORB,
  DISC1_RELN = A1$DISC1_RELN,
  STMN2 = C$STMN2
)

# Pivot longer for plotting
SP_stmn2_long <- SP_stmn2 %>%
  pivot_longer(cols = c(TLE4_CCBE1, CUX2_RORB, DISC1_RELN),
               names_to = "cell_type",
               values_to = "proportion")

# Calculate Spearman correlations
corr_stats_stmn2 <- SP_stmn2_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(STMN2, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(STMN2, proportion, method = "spearman")$p.value
  )

# Create labels
corr_stats_stmn2$label <- paste0("Rho = ", round(corr_stats_stmn2$rho, 2), 
                                 "\nP = ", signif(corr_stats_stmn2$pval, 2))

# Adjust text position manually
corr_stats_stmn2$x <- 0.13
corr_stats_stmn2$y <- c(0.000013, 0.000009, 0.000005)

# Set colors
cell_colors <- c(
  "TLE4_CCBE1" = "#F8766D",
  "CUX2_RORB" = "#00BFC4",
  "DISC1_RELN" = "orchid"
)

# Plot
ggplot(SP_stmn2_long, aes(x = STMN2, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = cell_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs STMN2 expression",
    x = "STMN2 expression",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats_stmn2,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 4, hjust = 0)
#######################################################################################

#C9 cell proportions Corr plots

### ACSL3

# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  OPC = A1$OPC,
  GFAP_neg = A1$GFAP.neg,
  RORB_ADGRL4 = A1$RORB_ADGRL4,
  RORB_LRRK1 = A1$RORB_LRRK1,
  ACSL3 = C$ACSL3
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(OPC, RORB_ADGRL4, RORB_LRRK1, GFAP_neg),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(ACSL3, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(ACSL3, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.31
corr_stats$y <- c(0.22, 0.17, 0.12, 0.07)  # one per cell type, spaced out



# Define inverted colors
cell_colors <- c(
  "OPC" = "#F8766D",
  "RORB_ADGRL4" = "#00BFC4",
  "RORB_LRRK1" = "orchid",
  "GFAP_neg" = "orange"
)

# Plot with annotations
ggplot(SP_long, aes(x = ACSL3, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = cell_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs pTDP-43 density",
    x = "ACSL3 density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0)



### FociAnti

# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  PCP4_NXPH2 = A1$PCP4_NXPH2,
  fociAnti = C$fociANTI
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(PCP4_NXPH2),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(fociAnti, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(fociAnti, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.3  # adjust as needed
corr_stats$y <- c(0.00001)  # one y per cell type for visual spacing

# Define inverted colors
inverted_colors <- c(
  "PCP4_NXPH2" = "#00BFC4"     # cyan
)

# Plot with annotations
ggplot(SP_long, aes(x = fociAnti, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs C9 cell proportions",
    x = "fociAnti density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0)


### FociSENSE

# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  PCP4_NXPH2 = A1$PCP4_NXPH2,
  PVALB_MYBPC1 = A1$PVALB_MYBPC1,
  Pericyte = A1$Pericyte,
  VIP_LAMA3 = A1$VIP_LAMA3,
  capillary = A1$Capillary,
  RORB_LRRK1 = A1$RORB_LRRK1,
  PVALB_CEMIP = A1$PVALB_CEMIP,
  fociSense = C$fociSENSE
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(  PCP4_NXPH2 ,PVALB_MYBPC1,Pericyte ,VIP_LAMA3, capillary, RORB_LRRK1, PVALB_CEMIP),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(fociSense, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(fociSense, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.3  # adjust as needed
corr_stats$y <- seq(0.01, 0.15, length.out = nrow(corr_stats))

# Define inverted colors
inverted_colors <- c(
  "PCP4_NXPH2" = "#00BFC4",     # cyan
  "PVALB_MYBPC1" = "red",
  "Pericyte" = "green",
  "VIP_LAMA3" = "orange",
  "capillary" = "pink",
  "RORB_LRRK1" = "purple",
  "PVALB_CEMIP" = "magenta"
)

# Plot with annotations
ggplot(SP_long, aes(x = fociSense, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs C9 cell proportions",
    x = "fociSense density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0)




### LNCRNA

# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  Oligo = A1$Oligo,
  THEMIS_NR4A2 = A1$THEMIS_NR4A2,
  RORB_LRRK1 = A1$RORB_LRRK1,
  lncRNA = C$lncRNA
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(   Oligo ,THEMIS_NR4A2,RORB_LRRK1),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(lncRNA, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(lncRNA, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.22  # adjust as needed
corr_stats$y <- seq(0.01, 0.15, length.out = nrow(corr_stats))

# Define inverted colors
inverted_colors <- c(
  "Oligo" = "#00BFC4",     # cyan
  "THEMIS_NR4A2" = "red",
  "RORB_LRRK1" = "purple"
)

# Plot with annotations
ggplot(SP_long, aes(x = lncRNA, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs C9 cell proportions",
    x = "lncRNA density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0)

## PolyGA

# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  Oligo = A1$Oligo,
  THEMIS_TMEM233 = A1$THEMIS_TMEM233,
  PolyGA = C$polyGA
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(THEMIS_TMEM233),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(PolyGA, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(PolyGA, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.000022  # adjust as needed
corr_stats$y <- seq(0.0001, 0.15, length.out = nrow(corr_stats))

# Define inverted colors
inverted_colors <- c("THEMIS_TMEM233" = "red")

# Plot with annotations
ggplot(SP_long, aes(x = PolyGA, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs C9 cell proportions",
    x = "PolyGA density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0)



## STMN2

# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  Oligo = A1$Oligo,
  THEMIS_TMEM233 = A1$THEMIS_TMEM233,
  TLE4_CCBE1 = A1$TLE4_CCBE1,
  PolyGP = C$polyGP
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(THEMIS_TMEM233, TLE4_CCBE1),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(PolyGP, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(PolyGP, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.000012  # adjust as needed
corr_stats$y <- seq(0.00006, 0.000015, length.out = nrow(corr_stats))

# Define inverted colors
inverted_colors <- c("THEMIS_TMEM233" = "red",
                     "TLE4_CCBE1" = "blue")

# Plot with annotations
ggplot(SP_long, aes(x = PolyGP, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs C9 cell proportions",
    x = "PolyGP density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0)


### STMN2

# Load data
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/C9/theta.state_cellstate.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_C9_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)
A1 <- A[C$X,]

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  GFAP_pos = A1$GFAP.pos,
  THEMIS_TMEM233 = A1$THEMIS_TMEM233,
  TLE4_SEMA3D = A1$TLE4_SEMA3D,
  RORB_LRRK1 = A1$RORB_LRRK1,
  STMN2 = C$STMN2
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c( GFAP_pos, THEMIS_TMEM233, TLE4_SEMA3D, RORB_LRRK1),
               names_to = "cell_type",
               values_to = "proportion")

# Compute Spearman correlations
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    rho = cor(STMN2, proportion, method = "spearman", use = "complete.obs"),
    pval = cor.test(STMN2, proportion, method = "spearman")$p.value
  )

# Format labels
corr_stats$label <- paste0("Rho = ", round(corr_stats$rho, 2), 
                           "\nP = ", signif(corr_stats$pval, 2))

# Set positions for annotation
corr_stats$x <- 0.03  # adjust as needed
corr_stats$y <- seq(0.1, 0.18, length.out = nrow(corr_stats))

# Define inverted colors
inverted_colors <- c(
  "GFAP_pos" = "#00BFC4",     # cyan
  "THEMIS_TMEM233" = "red",
  "TLE4_SEMA3D" = "green",
  "RORB_LRRK1" = "purple"
)

# Plot with annotations
ggplot(SP_long, aes(x = STMN2, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs C9 cell proportions",
    x = "STMN2 density",
    y = "Cell proportion",
    color = "Cell Type"
  ) +
  geom_text(data = corr_stats,
            aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE,
            size = 5, hjust = 0)











## OLD ##
###################################################################################################

#CELL PROPS


# V1

A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/theta.state_original.csv", row.names = 1)
C <-  read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")

C$X <- gsub("long", "X", C$X)

rownames(A)

A1 <- A[C$X,]

SP <- data.frame(GFAP_pos = A1$GFAP.pos,
                 RORB_LRRK1 = A1$RORB_LRRK1,
                 pTDP43 = C$TDP43b)
library("ggpubr")

ggscatter(
  data = SP,
  x = "pTDP43",
  y = c("GFAP_pos", "RORB_LRRK1"),
  add = "reg.line",      # Adds a regression line
  conf.int = TRUE,       # Adds confidence interval
  cor.coef = TRUE,       # Adds correlation coefficient
  cor.method = "spearman",# Method for correlation
  xlab = "TDP-43 Levels",
  ylab = "Cell proportion",
  title = "Cell proportion vs pTDP-43 density"
)
## V2
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  GFAP_pos = A1$GFAP.pos,
  RORB_LRRK1 = A1$RORB_LRRK1,
  pTDP43 = C$TDP43b
)

# Convert to long format for ggplot
SP_long <- SP %>%
  pivot_longer(cols = c(GFAP_pos, RORB_LRRK1),
               names_to = "cell_type",
               values_to = "proportion")

# Compute correlation statistics for each cell type
corr_stats <- SP_long %>%
  group_by(cell_type) %>%
  summarise(
    cor = cor(proportion, pTDP43, method = "spearman"),
    pval = cor.test(proportion, pTDP43, method = "spearman")$p.value,
    .groups = "drop"
  )

# Merge correlation stats back to the data for plotting
corr_labels <- corr_stats %>%
  mutate(
    label = paste0("r = ", round(cor, 2), "\np = ", signif(pval, 2)),
    x = max(SP$pTDP43) * 0.98,
    y = c(max(SP$GFAP_pos), max(SP$RORB_LRRK1)) * 1.05
  )

# Plot
ggplot(SP_long, aes(x = pTDP43, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text(data = corr_labels, aes(x = x, y = y, label = label, color = cell_type),
            inherit.aes = FALSE, hjust = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs pTDP-43 density",
    x = "TDP-43 Levels",
    y = "Cell proportion",
    color = "Cell Type"
  )
## V3

library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data
SP <- data.frame(
  GFAP_pos = A1$GFAP.pos,
  RORB_LRRK1 = A1$RORB_LRRK1,
  pTDP43 = C$TDP43b
)

# Convert to long format
SP_long <- SP %>%
  pivot_longer(cols = c(GFAP_pos, RORB_LRRK1),
               names_to = "cell_type",
               values_to = "proportion")

# Define inverted colors
inverted_colors <- c(
  "GFAP_pos" = "#00BFC4",     # Originally RORB_LRRK1 color (cyan)
  "RORB_LRRK1" = "#F8766D"    # Originally GFAP_pos color (red)
)

# Plot without correlation labels
ggplot(SP_long, aes(x = pTDP43, y = proportion, color = cell_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = inverted_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cell proportion vs pTDP-43 density",
    x = "pTDP-43 density",
    y = "Cell proportion",
    color = "Cell Type"
  )

########################################################################################################################################

# NPTX2

A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)
C <-  read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")

C$X <- gsub("long", "X", C$X)

rownames(A)

A1 <- A[C$X,]

SP <- data.frame(NPTX2 = A1$NPTX2,
                 pTDP43 = C$TDP43b)

corr <- cor(SP$NPTX2, SP$pTDP43, method = "spearman")
# test_result <- cor.test(SP$NPTX2, SP$pTDP43, method = "spearman")
library("ggpubr")

ggscatter(
  data = SP,
  x = "NPTX2",
  y = "pTDP43",
  add = "reg.line",      # Adds a regression line
  conf.int = TRUE,       # Adds confidence interval
  cor.coef = TRUE,       # Adds correlation coefficient
  cor.method = "spearman",# Method for correlation
  xlab = "NPTX2 Expression",
  ylab = "TDP-43 Levels",
  title = "NPTX2 vs TDP-43"
)
## V2

library(ggplot2)

# Llegir dades
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")

C$X <- gsub("long", "X", C$X)
A1 <- A[C$X, ]

# Preparar dades
SP <- data.frame(
  NPTX2 = A1$NPTX2,
  pTDP43 = C$TDP43b
)

# Definir el color manualment (mateix que RORB_LRRK1 anterior)
nptx2_color <- "#F8766D"

# Fer la gràfica
ggplot(SP, aes(x = pTDP43, y = NPTX2)) +
  geom_point(color = nptx2_color, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = nptx2_color) +
  theme_minimal(base_size = 14) +
  labs(
    title = "NPTX2 vs TDP-43",
    x = "pTDP-43 density",
    y = "NPTX2 Expression"
  )

#################################################################################################
#BOXPLOTS


### EXAMPLE URI

CELLS <- read.csv('/home/usuario/Documentos/AFTD/data/Elements/BlackTotalRNA/DESeq2/CellSubClusters_ALL_AAIC.csv')

CELLS$Cell_type_reordered <- factor(CELLS$Cell_type, levels=c("Ex8", "Ex4", "Mic2", "Ast2", "Ast3", "Oli5"))

png("/home/usuario/Documentos/AFTD/data/Elements/BlackTotalRNA/DESeq2/AAIC_3.pdf", width=20)

ggplot(CELLS,
       aes(x=Phenotype,
           y=Cell_proportion,
           fill=Phenotype)) +
  (theme_classic()) +   
  labs(x="", 
       y="Cell type proportion (%)")+  
  geom_boxplot(size=0.4, 
               outlier.alpha = 0.3, 
               alpha=0.5)+   
  facet_wrap(~`Cell_type_reordered`, 
             scales = "free", 
             nrow=1)+   
  theme_bw(base_size=12)+   
  theme(strip.text=element_text(face="bold"), 
        axis.title=element_text(size=12))+   
  scale_fill_manual(values=c("darkgoldenrod2", 
                             "deepskyblue4", 
                             "red", 
                             "black"))+  
  theme(panel.spacing = unit(5, "point"))+ 
  theme(legend.position="right")+ 
  scale_x_discrete(limits=rev(levels(CELLS$Phenotype)))

###### dades

# metadatas
SP_metadata <- readxl::read_xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx")
case_legend <- read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt", sep =  "\t", row.names = 1)
case_legend <- case_legend[case_legend$GROUP != "FTLD-TDP-C",]
case_legend$GROUP[case_legend$GROUP != "Control"] <- "TDP"

# RORB_LRRK1 (Sant pau)
sp1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/theta.state_original.csv", row.names = 1)

# RORB_LRRK1 (Validació)
n1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/theta.state_cellstate.csv", row.names = 1)

# NPTX2 (Sant pau)
sp2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)

# NPTX2 (Validació)
n2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)



####################### Sant pau LRRK1 prop #####################################
library(ggplot2)

sp1_ids <- gsub("X", "", rownames(sp1))
SP_metadata_filtered <- SP_metadata[SP_metadata$sample.ID %in% sp1_ids, ]
SP_metadata_filtered <- SP_metadata_filtered[match(sp1_ids, SP_metadata_filtered$sample.ID), ]
Pheno <- SP_metadata_filtered$group.ID
Cell_proportion <- sp1[rownames(sp1) %in% paste0("X", SP_metadata_filtered$sample.ID), "RORB_LRRK1"]
stopifnot(length(Pheno) == length(Cell_proportion))

CELLS <- data.frame(
  Phenotype = Pheno,
  Cell_proportion = Cell_proportion
)

colnames(CELLS)<- c("Phenotype", "Cell_proportion")
ggplot(CELLS,
       aes(x=Phenotype,
           y=Cell_proportion,
           fill=Phenotype)) +
  (theme_classic()) +   
  labs(x="", 
       y="Cell type proportion (%)")+  
  geom_boxplot(size=0.4, 
               outlier.alpha = 0.3, 
               alpha=0.5)+   
  #  facet_wrap(~`Cell_type_reordered`, 
  #             scales = "free", 
  #             nrow=1)+   
  theme_bw(base_size=12)+   
  theme(strip.text=element_text(face="bold"), 
        axis.title=element_text(size=12))+   
  scale_fill_manual(values=c("darkgoldenrod2", 
                             "deepskyblue4", 
                             "red", 
                             "black"))+  
  theme(panel.spacing = unit(5, "point"))+ 
  theme(legend.position="right")+ 
  scale_x_discrete(limits=rev(levels(CELLS$Phenotype)))


kruskal.test(Cell_proportion ~ Phenotype, data = CELLS)


###########################NPTX2 STPAU#############################

sp2_ids <- gsub("X", "", rownames(sp2))
SP_metadata_filtered <- SP_metadata[SP_metadata$sample.ID %in% sp2_ids, ]
SP_metadata_filtered <- SP_metadata_filtered[match(sp2_ids, SP_metadata_filtered$sample.ID), ]
Pheno <- SP_metadata_filtered$group.ID
Cell_proportion <- sp2[rownames(sp2) %in% paste0("X", SP_metadata_filtered$sample.ID), "NPTX2"]

stopifnot(length(Pheno) == length(Cell_proportion))

CELLS <- data.frame(
  Phenotype = Pheno,
  Cell_proportion = Cell_proportion
)

CELLS$sample.ID <- SP_metadata_filtered$sample.ID
colnames(CELLS)<- c("Phenotype", "Cell_proportion")

ggplot(CELLS,
       aes(x=Phenotype,
           y=Cell_proportion,
           fill=Phenotype)) +
  (theme_classic()) +   
  labs(x="", 
       y="NPTX2 expression")+  
  geom_boxplot(size=0.4, 
               outlier.alpha = 0.3, 
               alpha=0.5)+   
  #  facet_wrap(~`Cell_type_reordered`, 
  #             scales = "free", 
  #             nrow=1)+   
  theme_bw(base_size=12)+   
  theme(strip.text=element_text(face="bold"), 
        axis.title=element_text(size=12))+   
  scale_fill_manual(values=c("darkgoldenrod2", 
                             "deepskyblue4", 
                             "red", 
                             "black"))+  
  theme(panel.spacing = unit(5, "point"))+ 
  theme(legend.position="right")+ 
  scale_x_discrete(limits=rev(levels(CELLS$Phenotype)))


########################## RORB LRRK1 NEW##################################

n1_ids <- gsub("\\.", "-", rownames(n1))
common_ids <- intersect(rownames(case_legend), n1_ids)
case_legend_filtered <- case_legend[match(common_ids, rownames(case_legend)), ]
n1_filtered <- n1[match(common_ids, n1_ids), ]
Pheno <- case_legend_filtered$GROUP
Cell_proportion <- n1_filtered$RORB_LRRK1

CELLS <- data.frame(
  Phenotype = Pheno,
  Cell_proportion = Cell_proportion,
  sample.ID = common_ids  # optional, for checking alignment
)


ggplot(CELLS,
       aes(x=Phenotype,
           y=Cell_proportion,
           fill=Phenotype)) +
  (theme_classic()) +   
  labs(x="", 
       y="Cell type proportion (%)")+  
  geom_boxplot(size=0.4, 
               outlier.alpha = 0.3, 
               alpha=0.5)+   
  #  facet_wrap(~`Cell_type_reordered`, 
  #             scales = "free", 
  #             nrow=1)+   
  theme_bw(base_size=12)+   
  theme(strip.text=element_text(face="bold"), 
        axis.title=element_text(size=12))+   
  scale_fill_manual(values=c("darkgoldenrod2", 
                             "deepskyblue4", 
                             "red", 
                             "black"))+  
  theme(panel.spacing = unit(5, "point"))+ 
  theme(legend.position="right")+ 
  scale_x_discrete(limits=rev(levels(CELLS$Phenotype)))


########################## NPTX2 NEW##################################

n2_ids <- gsub("\\.", "-", rownames(n2))
common_ids <- intersect(rownames(case_legend), n2_ids)
case_legend_filtered <- case_legend[match(common_ids, rownames(case_legend)), ]
n2_filtered <- n2[match(common_ids, n2_ids), ]
Pheno <- case_legend_filtered$GROUP
Cell_proportion <- n2_filtered$NPTX2

CELLS <- data.frame(
  Phenotype = Pheno,
  Cell_proportion = Cell_proportion,
  sample.ID = common_ids
)

ggplot(CELLS,
       aes(x=Phenotype,
           y=Cell_proportion,
           fill=Phenotype)) +
  (theme_classic()) +   
  labs(x="", 
       y="NPTX2 expression")+  
  geom_boxplot(size=0.4, 
               outlier.alpha = 0.3, 
               alpha=0.5)+   
  #  facet_wrap(~`Cell_type_reordered`, 
  #             scales = "free", 
  #             nrow=1)+   
  theme_bw(base_size=12)+   
  theme(strip.text=element_text(face="bold"), 
        axis.title=element_text(size=12))+   
  scale_fill_manual(values=c("darkgoldenrod2", 
                             "deepskyblue4", 
                             "red", 
                             "black"))+  
  theme(panel.spacing = unit(5, "point"))+ 
  theme(legend.position="right")+ 
  scale_x_discrete(limits=rev(levels(CELLS$Phenotype)))

#pTDP-43

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# ----------------------
# 1. Read covariate data
# ----------------------
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)  # Ensure matching sample IDs

# -----------------------------------------------
# 2. Define cell types and genes you want to load
# -----------------------------------------------
cell_types <- c("RORB_POU3F2", "SMC", "RORB_LRRK1", "GFAP-neg", "Capillary", "T_Cell")

# Define folder
folder_path <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE ORIGINAL/"

# -------------------------------
# 3. Load and bind only needed files
# -------------------------------
expression_list <- lapply(cell_types, function(cell_type) {
  file_path <- paste0(folder_path, cell_type, ".csv")
  expr_df <- read.csv(file_path, row.names = 1)
  expr_df$ID <- rownames(expr_df)
  expr_long <- pivot_longer(expr_df, cols = -ID, names_to = "gene", values_to = "SP")
  expr_long$cell_type <- cell_type
  return(expr_long)
})

# Combine into one dataframe
expr_data <- bind_rows(expression_list)

# ----------------------------
# 4. Merge with pTDP-43 values
# ----------------------------
expr_data <- expr_data %>%
  left_join(C %>% select(X, TDP43b), by = c("ID" = "X")) %>%
  rename(pTDP43 = TDP43b)

# ----------------------------
# 5. Filter genes of interest
# ----------------------------
genes_of_interest <- c("KRT5", "CCDC184", "CHI3L1", "NRN1", "CST3", "TMEM160", "VGF", "NPTX2", "ADRA1D")
expr_data_filtered <- expr_data %>% filter(gene %in% genes_of_interest)

# -----------------------------
# 6. Plot using facet_wrap
# -----------------------------
ggplot(expr_data_filtered, aes(x = pTDP43, y = SP)) +
  geom_point(aes(color = gene), size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5) +
  facet_wrap(~ cell_type + gene, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Expression of Cell-State Markers vs. pTDP-43",
    x = "pTDP-43 density",
    y = "Gene expression (SP)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none"
  )

## STMN2

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# ----------------------
# 1. Read covariate data
# ----------------------
C <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/CORRELATION/COVARIABLES/FTD_TDP_neuropath_SOM.csv")
C$X <- gsub("long", "X", C$X)  # Ensure matching sample IDs

# -----------------------------------------------
# 2. Define cell types and genes you want to load
# -----------------------------------------------
cell_types <- c("RORB_POU3F2", "SMC", "RORB_LRRK1", "GFAP-neg", "Capillary", "T_Cell")

# Define folder
folder_path <- "/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE ORIGINAL/"

# -------------------------------
# 3. Load and bind only needed files
# -------------------------------
expression_list <- lapply(cell_types, function(cell_type) {
  file_path <- paste0(folder_path, cell_type, ".csv")
  expr_df <- read.csv(file_path, row.names = 1)
  expr_df$ID <- rownames(expr_df)
  expr_long <- pivot_longer(expr_df, cols = -ID, names_to = "gene", values_to = "SP")
  expr_long$cell_type <- cell_type
  return(expr_long)
})

# Combine into one dataframe
expr_data <- bind_rows(expression_list)

# ----------------------------
# 4. Merge with pTDP-43 values
# ----------------------------
expr_data <- expr_data %>%
  left_join(C %>% select(X, TDP43b), by = c("ID" = "X")) %>%
  rename(pTDP43 = TDP43b)

# ----------------------------
# 5. Filter genes of interest
# ----------------------------
genes_of_interest <- c(
  "KRT5", "CCDC184", "CHI3L1", "NRN1", "CST3", "TMEM160",
  "VGF", "NPTX2", "ADRA1D", "SOX11", "SLC24A4", "USH1C",
  "CPLX1", "C1QTNF4", "CCDC85B", "SCN9A", "RCAN2", "TAC1",
  "FOXJ1", "OR7D2", "FOXS1"
)

expr_data_filtered <- expr_data %>%
  filter(gene %in% genes_of_interest)

# -----------------------------
# 6. Plot using facet_wrap by cell type only
# -----------------------------
ggplot(expr_data_filtered, aes(x = pTDP43, y = SP, color = gene)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, size = 0.5, aes(group = gene), color = "black") +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Expression of Cell-State Markers vs. pTDP-43",
    x = "pTDP-43 density",
    y = "Gene expression (SP)",
    color = "Gene"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  )
###########################################################
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load data
counts <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/theta.state_original.csv", row.names = 1)
rownames(counts) <- gsub("X", "", rownames(counts))

meta <- read.xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx", 
                  row.names = 1, sheetIndex = 1)

# Filter metadata for groups of interest
meta <- meta[meta$group.ID %in% c("Healthy", "TDP"), ]
meta$group2.ID <- NULL

# Normalize counts rows to proportions summing to 1
counts_norm <- counts / rowSums(counts)

# Join counts and metadata by sample IDs
counts_norm$SampleID <- rownames(counts_norm)
meta$SampleID <- rownames(meta)

data_combined <- counts_norm %>%
  inner_join(meta, by = "SampleID") %>%
  filter(!is.na(group.ID))

# Sort samples by group.ID, keep order for plotting
data_combined <- data_combined %>%
  arrange(group.ID, SampleID) %>%
  mutate(SampleID = factor(SampleID, levels = SampleID))  # preserve order in plot

# Pivot longer for plotting
df_long <- data_combined %>%
  pivot_longer(cols = -c(SampleID, group.ID), 
               names_to = "CellState", 
               values_to = "Proportion")

# Define colors for groups
group_colors <- c("Healthy" = "blue", "TDP" = "red")

# Plot horizontal stacked barplot per sample
ggplot(df_long, aes(x = SampleID, y = Proportion, fill = CellState)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Sample", y = "Proportion of cells") +
  theme_minimal() +
  # Color y-axis labels by condition
  theme(axis.text.y = element_text(color = group_colors[df_long$group.ID][match(levels(df_long$SampleID), df_long$SampleID)])) +
  # Add manual legend for sample label colors
  guides(fill = guide_legend(title = "Cell State"))






library(dplyr)
library(tidyr)
library(ggplot2)

# Load data
counts <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/theta.state_cellstate.csv", row.names = 1)
rownames(counts) <- gsub("\\.", "-", rownames(counts))

meta <- read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt", row.names = 1)

# Filter metadata for groups of interest
meta <- meta[meta$GROUP != "FTLD-TDP-C", ]
meta$GROUP[meta$GROUP == "Control"] <- "Healthy"
meta$GROUP[meta$GROUP != "Healthy"] <- "TDP"
meta <- meta[meta$GROUP %in% c("Healthy", "TDP"), ]
meta <- meta[, "GROUP", drop = FALSE]

# Normalize counts rows to proportions summing to 1
counts_norm <- counts / rowSums(counts)

# Join counts and metadata by sample IDs
counts_norm$SampleID <- rownames(counts_norm)
meta$SampleID <- rownames(meta)

data_combined <- counts_norm %>%
  inner_join(meta, by = "SampleID") %>%
  filter(!is.na(GROUP)) %>%
  arrange(GROUP, SampleID) %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID)))  # preserve order

# Pivot longer for plotting
df_long <- data_combined %>%
  pivot_longer(cols = -c(SampleID, GROUP), names_to = "CellState", values_to = "Proportion")

# Define colors for groups
group_colors <- c("Healthy" = "blue", "TDP" = "red")

# Create named vector of colors for y-axis labels (sample names)
label_colors <- setNames(group_colors[data_combined$GROUP], data_combined$SampleID)

# Plot horizontal stacked barplot per sample
ggplot(df_long, aes(x = SampleID, y = Proportion, fill = CellState)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Sample", y = "Proportion of cells") +
  theme_minimal() +
  theme(axis.text.y = element_text(color = label_colors)) +
  guides(fill = guide_legend(title = "Cell State"))


#######################################################################
# Cell props New



CELLS <- data.frame(
  Phenotype = Pheno,
  Cell_proportion = Cell_proportion,
  sample.ID = common_ids
)

ggplot(CELLS,
       aes(x=Phenotype,
           y=Cell_proportion,
           fill=Phenotype)) +
  (theme_classic()) +   
  labs(x="", 
       y="NPTX2 expression")+  
  geom_boxplot(size=0.4, 
               outlier.alpha = 0.3, 
               alpha=0.5)+   
  #  facet_wrap(~`Cell_type_reordered`, 
  #             scales = "free", 
  #             nrow=1)+   
  theme_bw(base_size=12)+   
  theme(strip.text=element_text(face="bold"), 
        axis.title=element_text(size=12))+   
  scale_fill_manual(values=c("darkgoldenrod2", 
                             "deepskyblue4", 
                             "red", 
                             "black"))+  
  theme(panel.spacing = unit(5, "point"))+ 
  theme(legend.position="right")+ 
  scale_x_discrete(limits=rev(levels(CELLS$Phenotype)))


SP_metadata <- readxl::read_xlsx("/media/jaumatell/datos/URI/BAYESPRISM_12_3/FTLD_BULK/METADATA/decoder_DeSeq2_FTD_FINAL.xlsx")
case_legend <- read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt", sep =  "\t", row.names = 1)
case_legend <- case_legend[case_legend$GROUP != "FTLD-TDP-C",]
case_legend$GROUP[case_legend$GROUP != "Control"] <- "TDP"

# RORB_LRRK1 (Sant pau)
sp1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/theta.state_original.csv", row.names = 1)

# RORB_LRRK1 (Validació)
n1<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/theta.state_cellstate.csv", row.names = 1)

# NPTX2 (Sant pau)
sp2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/FTLD/TDP/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)

# NPTX2 (Validació)
n2<-read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)




case_legend <- read.delim("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt", sep = "\t", row.names = 1)
case_legend <- case_legend[case_legend$GROUP != "FTLD-TDP-C", ]
case_legend$GROUP[case_legend$GROUP != "Control"] <- "TDP"

n1 <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/theta.state_cellstate.csv", row.names = 1)
n2 <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/BAYESPRISM/NEW/CELL STATE ORIGINAL/RORB_LRRK1.csv", row.names = 1)

n2_ids <- gsub("\\.", "-", rownames(n2))
common_ids <- intersect(rownames(case_legend), n2_ids)
case_legend_filtered <- case_legend[match(common_ids, rownames(case_legend)), ]
n2_filtered <- n2[match(common_ids, n2_ids), ]
Pheno <- case_legend_filtered$GROUP

target_cell_states <- c(
  "VAT1L_EYA4", "PVALB_CEMIP", "Arterial", "CUX2_RORB", "PCP4_NXPH2",
  "GFAP.pos", "SST_ADAMTS19", "PVALB_MYBPC1", "CDH4_CCK", "LAMP5_PMEPA1",
  "PVALB_PTHLH", "DISC1_CCK", "RORB_POU3F2", "SST_GALNT14", "THEMIS_TMEM233",
  "RORB_LRRK1", "DISC1_RELN", "LAMP5_CA3", "VIP_HTR2C", "VIP_LAMA3", "SST_BRINP3"
)

nice_names <- c(
  "VAT1L EYA4", "PVALB CEMIP", "Arterial", "CUX2 RORB", "PCP4 NXPH2",
  "GFAP+", "SST ADAMTS19", "PVALB MYBPC1", "CDH4 CCK", "LAMP5 PMEPA1",
  "PVALB PTHLH", "DISC1 CCK", "RORB POU3F2", "SST GALNT14", "THEMIS TMEM233",
  "RORB LRRK1", "DISC1 RELN", "LAMP5 CA3", "VIP HTR2C", "VIP LAMA3", "SST BRINP3"
)

cell_data <- data.frame()

for (i in seq_along(target_cell_states)) {
  cell_type <- target_cell_states[i]
  label <- nice_names[i]
  cell_prop <- n1[rownames(n1) %in% gsub("-" , "\\.", rownames(case_legend_filtered)), cell_type]
  temp_df <- data.frame(
    Phenotype = Pheno,
    Cell_proportion = cell_prop,
    Cell_state = label
  )
  cell_data <- rbind(cell_data, temp_df)
}

cell_data$Phenotype <- factor(cell_data$Phenotype, levels = c("Control", "TDP"))

library(dplyr)
library(ggplot2)

max_y_per_state <- cell_data %>%
  group_by(Cell_state) %>%
  summarise(y = max(Cell_proportion, na.rm = TRUE) * 0.95)

annotation_df <- data.frame(
  Cell_state = c(
    "VAT1L EYA4", "PVALB CEMIP", "Arterial", "CUX2 RORB", "PCP4 NXPH2",
    "GFAP+", "SST ADAMTS19", "PVALB MYBPC1", "CDH4 CCK", "LAMP5 PMEPA1",
    "PVALB PTHLH", "DISC1 CCK", "RORB POU3F2", "SST GALNT14", "THEMIS TMEM233",
    "RORB LRRK1", "DISC1 RELN", "LAMP5 CA3", "VIP HTR2C", "VIP LAMA3", "SST BRINP3"
  ),
  label = c(
    "P = 0.0001\nFC = 0.04", "P = 0.0004\nFC = 0.60", "P = 0.0007\nFC = 1.65", "P = 0.0012\nFC = 0.50", "P = 0.0017\nFC = 0.58",
    "P = 0.0032\nFC = 115.25", "P = 0.0046\nFC = 0.71", "P = 0.0046\nFC = 0.62", "P = 0.0054\nFC = 0.60", "P = 0.0082\nFC = 0.97",
    "P = 0.0085\nFC = 0.69", "P = 0.0092\nFC = 0.67", "P = 0.0099\nFC = 0.62", "P = 0.0119\nFC = 0.70", "P = 0.0137\nFC = 0.65",
    "P = 0.0164\nFC = 0.63", "P = 0.0239\nFC = 0.72", "P = 0.0256\nFC = 0.62", "P = 0.0291\nFC = 0.67", "P = 0.0331\nFC = 0.64", "P = 0.0353\nFC = 0.73"
  ),
  x = 1.5
)

annotation_df <- left_join(annotation_df, max_y_per_state, by = "Cell_state")

ggplot(cell_data, aes(x = Phenotype, y = Cell_proportion, fill = Phenotype)) +
  geom_boxplot(size = 0.4, outlier.alpha = 0.3, alpha = 0.5) +
  facet_wrap(~Cell_state, scales = "free", nrow = 3) +
  labs(x = "", y = "Cell type proportion (%)") +
  scale_fill_manual(values = c("Control" = "deepskyblue4", "TDP" = "darkgoldenrod2")) +
  geom_text(data = annotation_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.5) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    panel.spacing = unit(5, "pt"),
    legend.position = "right"
  )

library("ggplot2")
library("openxlsx")
library("dplyr")
library("coin")
library("AnnotationDbi")
library("org.Hs.eg.db")
library('EnsDb.Hsapiens.v86')
library("Seurat")
library("edgeR")
library("GO.db")
library("xlsx")
suppressWarnings(library("BayesPrism"))
suppressWarnings(library(org.Hs.eg.db))
library(dplyr) 
library("AnnotationDbi")
library("org.Hs.eg.db")
library('EnsDb.Hsapiens.v86')
library("Seurat")

################################################################################
# Single cell data
load("/media/jaumatell/datos/URI/BAYESPRISM_12_3/SINGLE_CELL/DATA/Merge_FC_complete.RData")
# Bulk data
bk.dat <-  read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/DATA/merged_gene_count_FCX.csv")
###############################################################################
#Set colnames and rownames
bk.dat<-bk.dat[, -c(1,2,3,4,5,7)]
bk.dat <- aggregate(. ~ GeneName, data = bk.dat, FUN = sum) # Aggregate rows with same gene
rownames(bk.dat)<-bk.dat[,1]
bk.dat <- bk.dat[,-1]
bk.dat <- t(bk.dat)


sc.dat <- t(merged@assays$Assay_name$counts)
# El matrix es el que dona problemes!!!!!!
#sc.dat <- t(as.matrix(merged@assays$Assay_name$counts))
A <- read.csv("/media/jaumatell/datos/URI/BAYESPRISM_12_3/SINGLE_CELL/METADATA/S03_annotation_mapping.csv")
A[44,"SubType" ] <- "T_Cell"
A[20,"SubType" ] <- "OPC"
A[19,"SubType" ] <- "Oligo"
A[18,"SubType" ] <- "Micro"


# Update metadata to take into account DGE groups.

merged@meta.data <- merged@meta.data %>%
  left_join(A, by = c("cellstate" = "SubType"))


# Cell type and state
cell.state.labels <- merged@meta.data$cellstate
cell.type.labels <- merged@meta.data$CellType
rm(merged)

# Pre-processing
plot.cor.phi (input=sc.dat, 
              input.labels=as.factor(cell.state.labels),
              title="cell state correlation",
              pdf.prefix="gbm.cor.cs",
              cexRow=0.2, 
              cexCol=0.2,
              margins=c(2,2)
)

plot.cor.phi (input=sc.dat,
              input.labels=as.factor(cell.type.labels),
              title="cell type correlation",
              pdf.prefix="gbm.cor.ct",
              cexRow=0.5, 
              cexCol=0.5,
)

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=as.factor(cell.type.labels),
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE, #return the data used for plotting.
  pdf.prefix="SC_stats"
)
write.csv(sc.stat, "sc.stat.csv")

print("bk.stat")
bk.stat <- plot.bulk.outlier(
  bulk.input = bk.dat,
  sc.input = sc.dat,
  cell.type.labels = cell.type.labels,
  species = "hs",
  return.raw = TRUE,
  pdf.prefix = "BK_stats"
)
write.csv(bk.dat, "bk.stat.csv")


print("sc.stat.filtered")
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5
)
write.csv(sc.dat.filtered, "sc.dat.filtered.csv")

print("bk.stat.filtered")
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat,
                 pdf.prefix="Bulk_vs_Sc"
)

print("sc.dat.filtered.pc")
sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")
write.csv(sc.dat.filtered.pc, "sc.dat.filtered.pc.csv")

# SORT bk.dat so the colnames are sorted equaly in both tables. 
common_columns <- intersect(colnames(bk.dat), colnames(sc.dat.filtered.pc))
bk.dat <- bk.dat[, common_columns]
sc.dat.filtered.pc <- sc.dat.filtered.pc[, common_columns]

rm(sc.dat)
rm(sc.dat.filtered)
rm(bk.stat)
rm(sc.stat)

myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=as.matrix(bk.dat[,2:ncol(bk.dat)]),
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
save(myPrism, file = "myPrism.RData")

bp.res <- run.prism(prism = myPrism, n.cores=50,)
save(bp.res, file = "bp.res.RData")

theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
write.csv(theta, "theta.csv")

theta.cv <- bp.res@posterior.theta_f@theta.cv
write.csv(theta.cv, "theta.cv.csv")

bp.res.initial <- run.prism(prism = myPrism, 
                            n.cores=50, 
                            update.gibbs=FALSE)
save(bp.res.initial, file = "bp.res.initial.RData")
bp.res.update <- update.theta (bp = bp.res.initial)
save(bp.res.update, file = "bp.res.update.RData")

theta.state <- get.fraction (bp=bp.res.initial,
                             which.theta="first",
                             state.or.type="state")
write.csv(theta.state, "theta.state_cellstate.csv")

# Make new directories
new_dir_name <- "CELL TYPE"
if (file.exists(new_dir_name)){
  print("", end = "")
} else {
  dir.create(new_dir_name)
}

new_dir_name <- "CELL STATE ORIGINAL"
if (file.exists(new_dir_name)){
  print("", end = "")
} else {
  dir.create(new_dir_name)
}

for (cell_type in levels(as.factor(cell.type.labels))){
  exp_type <- get.exp(bp = bp.res, 
                      state.or.type = "type", 
                      cell.name = cell_type)
  write.csv(exp_type, paste0("CELL TYPE/", cell_type, ".csv"))
}


for (cell_state in levels(as.factor(cell.state.labels))){
  exp_type <- get.exp(bp = bp.res, 
                      state.or.type = "state", 
                      cell.name = cell_state)
  write.csv(exp_type, paste0("CELL STATE ORIGINAL/", cell_state, ".csv"))
}


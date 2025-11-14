library("sva")
library("DESeq2")

coldataTAU <- read.delim('/media/jaumatell/datos/URI/PROJECTE_SEURAT_BP/RiMod/rimod_ftd_dataset_table_v3.txt', sep = "\t", row.names = 1)
ctsTAU <- as.matrix(read.csv('/media/jaumatell/datos/URI/PROJECTE_SEURAT_BP/RiMod/rnaseq_salmon_results/counts.csv',sep=",", row.names = 1))

common_IDs <- intersect(rownames(coldataTAU), colnames(ctsTAU))

coldataTAU <- coldataTAU[common_IDs,]
ctsTAU <- ctsTAU[, common_IDs]

# PROCESS 
ddsTAU <- DESeqDataSetFromMatrix(countData = round(ctsTAU), colData = coldataTAU, design = ~ DiseaseCode)

idx <- rowSums( counts(ddsTAU, normalized=FALSE) >= 10 ) >= 20

ddsTAU_Filter <- ddsTAU[idx,]

ddsTAU_Filter_ <- DESeq(ddsTAU_Filter)

dat <- counts(ddsTAU_Filter_, normalized=TRUE)

mod <- model.matrix(~ DiseaseCode, colData(ddsTAU_Filter_))

mod0 <- model.matrix(~ 1, colData(ddsTAU_Filter_))

num.sv(dat, mod) #Aquí ens diu quants cal aplicar#

svseqTau <- svaseq(dat, mod, mod0)

ddsTAU_Filter_SV2 <- ddsTAU_Filter_

ddsTAU_Filter_SV2$SV1 <- svseqTau$sv[,1]

#ddsTAU_Filter_SV2$SV2 <- svseqTau$sv[,2]


###ddsTAU_Filter_SV2 conté ID, SV1 i SV2 (els dos SV que necessitem) ###

sva_dataframe_IDs <- rownames(coldataTAU)
sva_dataframe_disease <- coldataTAU[,"DiseaseCode"]
sva_dataframe_sv1 <- svseqTau$sv[,1]

sva_dataframe <- data.frame(
  ID = sva_dataframe_IDs,
  Disease = sva_dataframe_disease,
  SV1 = sva_dataframe_sv1,
  stringsAsFactors = FALSE
)

write.table(sva_dataframe, file = "/media/jaumatell/datos/URI/BAYESPRISM_12_3/RIMOD_BULK/METADATA/SVA.csv", sep = "\t", row.names = FALSE, quote = FALSE)

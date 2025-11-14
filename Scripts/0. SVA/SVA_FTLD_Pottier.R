library("sva")
library("DESeq2")

coldataTAU <- read.delim('/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/Sample_info.txt', sep = "\t", row.names = 1)
ctsTAU <- as.matrix(read.csv('/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/DATA/merged_gene_count_FCX.csv',sep=",", row.names = 1))
ctsTAU <- ctsTAU[,-c(1,2,3,4,5,6)]

coldataTAU <- coldataTAU[coldataTAU$GROUP!= "FTLD-TDP-C", ]

colnames(ctsTAU) <- gsub("\\.", "-", colnames(ctsTAU))
common_IDs <- intersect(rownames(coldataTAU), colnames(ctsTAU))

coldataTAU <- coldataTAU[common_IDs,]
ctsTAU <- ctsTAU[, common_IDs]

ctsTAU <- apply(ctsTAU, 2, as.numeric)


# PROCESS 
ddsTAU <- DESeqDataSetFromMatrix(countData = round(ctsTAU), colData = coldataTAU, design = ~ GROUP)

idx <- rowSums( counts(ddsTAU, normalized=FALSE) >= 10 ) >= 20

ddsTAU_Filter <- ddsTAU[idx,]

ddsTAU_Filter_ <- DESeq(ddsTAU_Filter)

dat <- counts(ddsTAU_Filter_, normalized=TRUE)

mod <- model.matrix(~ GROUP, colData(ddsTAU_Filter_))

mod0 <- model.matrix(~ 1, colData(ddsTAU_Filter_))

num.sv(dat, mod) #Aquí ens diu quants cal aplicar#

svseqTau <- svaseq(dat, mod, mod0)

ddsTAU_Filter_SV2 <- ddsTAU_Filter_

ddsTAU_Filter_SV2$SV1 <- svseqTau$sv[,1]

sva_dataframe_IDs <- rownames(coldataTAU)
sva_dataframe_disease <- coldataTAU[,"GROUP"]
sva_dataframe_sv1 <- svseqTau$sv[,1]

sva_dataframe <- data.frame(
  ID = sva_dataframe_IDs,
  Disease = sva_dataframe_disease,
  SV1 = sva_dataframe_sv1,
  stringsAsFactors = FALSE
)

write.table(sva_dataframe, file = "/media/jaumatell/datos/URI/BAYESPRISM_12_3/NEW_BULK/METADATA/SVA_correcte.csv", sep = "\t", row.names = FALSE, quote = FALSE)
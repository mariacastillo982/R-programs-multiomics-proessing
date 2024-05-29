
library("readxl")
library(dplyr)
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library("GEOquery")
library(genefilter)
library(clusterProfiler)
library(clariomdhumanhsrefseq.db)
library(clariomdhumanhsrefseqcdf)
library(pd.clariomdhuman.hs.refseq)
#transcriptomic_data=read.csv('/Users/mariacastillo/Desktop/HEATSTROKE/DATA CEL FILES/ALL/filtered_normalized_exprs_Symbol_timepoint_Groupname.csv')
transcriptomic_data= data.frame(exprs)
transcriptomic_data <- transcriptomic_data[, !grepl("CD", names(transcriptomic_data))]
#rownames(transcriptomic_data) = transcriptomic_data$X
#transcriptomic_data$X <- NULL

proteomics_data=read.csv('/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/PBMC_1st_Injection_DIANN_results20072023/normalized_data_Standard_Deviation.csv')
data_ = AnnotationDbi::select(org.Hs.eg.db, keys =proteomics_data$X, columns = "SYMBOL", keytype = "UNIPROT")
colnames(data_)[1] <- "X"
proteomics_data2 <- left_join(data_,proteomics_data)
#rownames(proteomics_data) = proteomics_data$X
#proteomics_data$X <- NULL

common_rows <- intersect(proteomics_data2$SYMBOL, rownames(transcriptomic_data))
indexes <- which(proteomics_data2$SYMBOL %in% common_rows, arr.ind = TRUE)
trans <- transcriptomic_data[common_rows, ]
trans$HS.02.Stroke.T1 = NULL
prot <- proteomics_data2[indexes, ]
duplicates <- prot$SYMBOL[duplicated(prot$SYMBOL)]
prot <- prot[!duplicated(prot$SYMBOL), ]
rownames(prot) = prot$SYMBOL
prot$X <- NULL
prot$SYMBOL = NULL
column_order_prot <- names(prot)
trans <- trans[column_order_prot]

cor.test(prot$HS.01.Stroke.T0,trans$HS.01.Stroke.T0)
cor.test(prot$HS.01.Stroke.T1,trans$HS.01.Stroke.T1)
cor.test(prot$HS.01.Stroke.T0,trans$HS.01.Stroke.T1)
cor.test(prot$HS.01.Stroke.T1,trans$HS.01.Stroke.T0)

cor.test(prot$HS.03.Stroke.T0,trans$HS.03.Stroke.T0)
cor.test(prot$HS.03.Stroke.T1,trans$HS.03.Stroke.T1)
cor.test(prot$HS.03.Stroke.T0,trans$HS.03.Stroke.T1)
cor.test(prot$HS.03.Stroke.T1,trans$HS.03.Stroke.T0)

cor.test(prot$HS.07.Stroke.T0,trans$HS.07.Stroke.T0)
cor.test(prot$HS.07.Stroke.T1,trans$HS.07.Stroke.T1)
cor.test(prot$HS.07.Stroke.T0,trans$HS.07.Stroke.T1)
cor.test(prot$HS.07.Stroke.T1,trans$HS.07.Stroke.T0)

cor.test(prot$HS.19.Stroke.T0,trans$HS.19.Stroke.T0)
cor.test(prot$HS.19.Stroke.T1,trans$HS.19.Stroke.T1)
cor.test(prot$HS.19.Stroke.T0,trans$HS.19.Stroke.T1)
cor.test(prot$HS.19.Stroke.T1,trans$HS.19.Stroke.T0)
write.csv(correlation, "/Users/mariacastillo/Desktop/correlation_all_genes.csv")

trans_t0=trans[, grepl("T0", names(trans))]
trans_t1=trans[, grepl("T1", names(trans))]

prot_t0=prot[, grepl("T0", names(prot))]
prot_t1=prot[, grepl("T1", names(prot))]

correlation_t0_t0 = cor(t(prot_t0),t(trans_t0))
correlation_t1_t1 = cor(t(prot_t1),t(trans_t1))
trans_t0$HS.14.Stroke.T0 = NULL
trans_t0$HS.02.Stroke.T0 = NULL
trans_t0$HS.16.Stroke.T0 = NULL
prot_t0$HS.02.Stroke.T0 = NULL
prot_t0$HS.14.Stroke.T0 = NULL
prot_t0$HS.16.Stroke.T0 = NULL
correlation_t0_t1 <- cor(t(prot_t0), t(trans_t1))


# Housekeeping genes
genes_int <- c("HPRT1", "RPLP0","UBC", "RPL13A",
                        "HSPA4L", "HSPB1", "HSP90AA1", "HSP90AB1", "HSPA1A", 
                      "HSPA1B", "HSPD1", "HSPA8", "HSP90B1", "HSPA6", 
                      "HSPA4", "HSPA9", "HSPA13", "HSPA2", "HSPE1", "HSPA14", 
                      "HSPH1", "DNAJB1","DNAJC1","DNAJC2", "DNAJC7", "DNAJC5", "DNAJA1", "DNAJA4")
trans_ss <- trans[genes_int, ]
prot_ss <- prot[genes_int, ]

correlation = cor(prot_ss,trans_ss)
write.csv(correlation, "/Users/mariacastillo/Desktop/correlation_target_genes.csv")

trans_t0=trans_ss[, grepl("T0", names(trans_ss))]
trans_t1=trans_ss[, grepl("T1", names(trans_ss))]
prot_t0=prot_ss[, grepl("T0", names(prot_ss))]
prot_t1=prot_ss[, grepl("T1", names(prot_ss))]

correlation_t0_t0 = cor(t(prot_t0),t(trans_t0))
correlation_t1_t1 = cor(t(prot_t1),t(trans_t1))
trans_t0$HS.14.Stroke.T0 = NULL
trans_t0$HS.02.Stroke.T0 = NULL
trans_t0$HS.16.Stroke.T0 = NULL
prot_t0$HS.02.Stroke.T0 = NULL
prot_t0$HS.14.Stroke.T0 = NULL
prot_t0$HS.16.Stroke.T0 = NULL
correlation_t1_t0 <- cor(t(trans_t1), t(prot_t0))
correlation_t0_t1 <- cor(t(trans_t0), t(prot_t1))

write.csv(correlation_t0_t0, "/Users/mariacastillo/Desktop/correlation_t0_t0.csv")
write.csv(correlation_t1_t1, "/Users/mariacastillo/Desktop/correlation_t1_t1.csv")
write.csv(correlation_t0_t1, "/Users/mariacastillo/Desktop/correlation_t0_t1.csv")
write.csv(correlation_t1_t0, "/Users/mariacastillo/Desktop/correlation_t1_t0.csv")

proteomic_DE=read.csv("/Users/mariacastillo/Desktop/HEATSTROKE/correlation transcriptome proteome/proteomic_DE.csv")
transcriptomic_DE=read.csv("/Users/mariacastillo/Desktop/HEATSTROKE/correlation transcriptome proteome/transcriptomic_DE.csv")
comp=intersect(proteomic_DE$SYMBOL, transcriptomic_DE$SYMBOL)
indexes <- which(proteomic_DE$SYMBOL %in% comp, arr.ind = TRUE)
proteomic_DE <- proteomic_DE[indexes, ]
proteomic_DE <- proteomic_DE[!duplicated(proteomic_DE$SYMBOL), ]
rownames(proteomic_DE)=proteomic_DE$SYMBOL
proteomic_DE$SYMBOL=NULL
indexes <- which(transcriptomic_DE$SYMBOL %in% comp, arr.ind = TRUE)
transcriptomic_DE <- transcriptomic_DE[indexes, ]
transcriptomic_DE <- transcriptomic_DE[!duplicated(transcriptomic_DE$SYMBOL), ]
rownames(transcriptomic_DE)=transcriptomic_DE$SYMBOL
transcriptomic_DE$SYMBOL=NULL
correlation_DE = cor(transcriptomic_DE,proteomic_DE)




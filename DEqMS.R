## ----LoadpackDisease--------------------------------------------------------------
library(DEqMS)
library(affy)
library(affycoretools)
library("readxl")
library("DEP")
library("dplyr")
library(org.Hs.eg.db)
library(readr)
library(stringr) 
library(dplyr)
library(M3C)
library(EnhancedVolcano)
library(clusterProfiler)
library('FactoMineR')
library(openxlsx)
library(ggforce)
library(ggplot2)
library(dplyr)
library(tidyr)

# PROTEINS
setwd("/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/PBMC_1st_Injection_DIANN_results20072023/")
meta <- read_excel("metadata.xlsx")
data<-readr::read_tsv("report.pg_matrix.tsv")
dat=data[6:59]
rownames(dat)=data$Protein.Group
colnames(dat)=meta$`Sample ID`

meta_Old <- meta[ meta$Disease=="Old", ]
meta_Young <- meta[ meta$Disease=="Young", ]
dat_Old = dat[,meta_Old$label]
rownames(dat_Old)=data$Protein.Group
dat_Young = dat[,meta_Young$label]
rownames(dat_Young)=data$Protein.Group

## ----log2transform1-----------------------------------------------------------
dat.log = log2(dat)
# Create a logical vector indicating which rows have NAs
na_rows = rowSums(is.na(dat.log)) > 0
#remove rows with NAs
dat.log = na.omit(dat.log)
protein_group_filtered = data$Protein.Group[!na_rows]
rownames(dat.log)=protein_group_filtered
## ----boxplot1-----------------------------------------------------------------
boxplot(dat.log,las=2,main="Cellular proteomic data: Proteins")
dev.copy(jpeg,filename="Barplot_before_Normalization_Proteins.jpg");
dev.off ();

# Identify the index of the reference protein in the data frame
reference_index <- which(rownames(dat.log) == 'P04406')
# Normalize each row by dividing it by the values of the reference row
vec=dat.log[reference_index,]
normalized_data=mapply('/', dat.log, vec)
rownames(normalized_data)=rownames(dat.log)
normalized_data=data.frame(normalized_data)
write.csv(normalized_data, "normalized_data_HKG.csv")

boxplot(normalized_data,las=2,main="Cellular proteomic data: Proteins")
dev.copy(jpeg,filename="Barplot_before_Normalization_Proteins_HKG.jpg");
dev.off ();

# Here the data is already median centered, we skip the following step. 
dat.log = equalMedianNormalization(dat.log)
boxplot(dat.log,las=2,main="Cellular proteomic data: Proteins")
dev.copy(jpeg,filename="Barplot_after_Normalization_Proteins.jpg");
dev.off ();
write.csv(dat.log, "normalized_data_Standard_Deviation.csv")

value=as.numeric(unlist(dat.log))
qqnorm(value,pch = 1, frame = FALSE)
qqline(value, col = "steelblue", lwd = 2)
dev.copy(jpeg,filename="Q-Q plot Proteins.jpg");
dev.off ();

## -------------PCA-------------------------------------------------------------
gs= factor(meta$condition)
Disease = factor(meta$Disease)
Gender = factor(meta$Gender)
PCA_raw <- prcomp(t(dat.log), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Disease = Disease,
                     Gender = Gender,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Disease) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Disease) +  # Color by Disease
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = TRUE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "PCA_Proteins_PC1_PC2_Disease.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Disease) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Disease) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = TRUE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "PCA_Proteins_PC2_PC3_Disease.jpg")
dev.off()

## ----design-------------------------------------------------------------------
x=meta$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)
# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
colnames(design) <- c("Stress","T0","T1")

fit1 <- lmFit(dat.log, design)
contrast =  makeContrasts(contrasts="Stress-T1",levels=design)
fit2 <- contrasts.fit(fit1, contrast)
fit3 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit3)
summary(DE_genes2)
# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit3, adjust="BH", sort.by ="p", number=Inf, p.value=0.05)
top_genes=na.omit(top_genes)
head(top_genes)

output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/TopGenes_DEP_norm.xlsx'
wb <- createWorkbook()

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "SYMBOL", keytype = "UNIPROT")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='UNIPROT'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, by = "UNIPROT", all = TRUE)
gene_list=na.omit(gene_list)

write.csv(gene_list, "DE_Proteins_Stress_T1_norm2.csv")
#output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/TopGenes_DEP.xlsx'

addWorksheet(wb, sheetName = 'Stress-T1 (Proteins)')
writeData(wb, sheet = 'Stress-T1 (Proteins)', x = gene_list)

EnhancedVolcano(top_genes,
                lab = rownames(top_genes),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = c(min(top_genes[['logFC']], na.rm = TRUE) - 0.5, max(top_genes[['logFC']], na.rm = TRUE) +0.5),
                ylim = c(0, max(-log10(top_genes[['adj.P.Val']]), na.rm = TRUE) + 5),
                title = 'Proteins: Control vs. Heat stroke T1',
                pCutoff = 0.0001,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="Volcano_Proteins_Stress_T1_norm.jpg");
dev.off ();

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "UNIPROT")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='UNIPROT'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, BH = "UNIPROT", all = TRUE)
gene_list=na.omit(gene_list)

# GENE ENRICHMENT ANALYSIS
OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_df <- as.data.frame(ego)
head(ego_df[1:7])
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 
dev.copy(jpeg,filename="Barplot_GEA_Proteins_Stress_T1_norm.jpg");
dev.off ();

contrast =  makeContrasts(contrasts="Stress-T0",levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)
# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit3, adjust="BH", sort.by ="p", number=Inf, p.value=0.05)
top_genes=na.omit(top_genes)

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "SYMBOL", keytype = "UNIPROT")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='UNIPROT'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, by = "UNIPROT", all = TRUE)
gene_list=na.omit(gene_list)

write.csv(gene_list, "DE_Proteins_Stress_T0_norm.csv")
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Stress-T0 (Proteins)')
writeData(wb, sheet = 'Stress-T0 (Proteins)', x = gene_list)


EnhancedVolcano(top_genes,
                lab = rownames(top_genes),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = c(min(top_genes[['logFC']], na.rm = TRUE) - 0.5, max(top_genes[['logFC']], na.rm = TRUE) +0.5),
                ylim = c(0, max(-log10(top_genes[['adj.P.Val']]), na.rm = TRUE) + 5),
                title = 'Proteins: Control vs. Heat stroke T0',
                pCutoff = 0.0001,
                FCcutoff = 1,
                pointSize = 1.0)

dev.copy(jpeg,filename="Volcano_Proteins_Stress_T0_norm.jpg");
dev.off ();

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "UNIPROT")
#top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='UNIPROT'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, BH = "UNIPROT", all = TRUE)
gene_list=na.omit(gene_list)

# GENE ENRICHMENT ANALYSIS
OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ego_df <- as.data.frame(ego)
head(ego_df[1:7])
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 
dev.copy(jpeg,filename="Barplot_GEA_Proteins_Stress_T0_norm.jpg");
dev.off ();

contrast =  makeContrasts(contrasts="T0-T1",levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)
# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit3, adjust="BH", sort.by ="p", number=Inf, p.value=0.05)
top_genes=na.omit(top_genes)
gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "SYMBOL", keytype = "UNIPROT")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='UNIPROT'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, by = "UNIPROT", all = TRUE)
gene_list=na.omit(gene_list)

write.csv(gene_list, "DE_Proteins_T0_T1_norm.csv")
addWorksheet(wb, sheetName = 'T0-T1 (Proteins)')
writeData(wb, sheet = 'T0-T1 (Proteins)', x = gene_list)

# VOLCANO PLOT
EnhancedVolcano(top_genes,
                lab = rownames(top_genes),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = c(min(top_genes[['logFC']], na.rm = TRUE) - 0.5, max(top_genes[['logFC']], na.rm = TRUE) +0.5),
                ylim = c(0, max(-log10(top_genes[['adj.P.Val']]), na.rm = TRUE) + 5),
                title = 'Proteins: Heat stroke T0 vs. Heat stroke T1',
                pCutoff = 0.0001,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="Volcano_Proteins_T0_T1_norm.jpg");
dev.off ();

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "UNIPROT")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='UNIPROT'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, BH = "UNIPROT", all = TRUE)
gene_list=na.omit(gene_list)

# GENE ENRICHMENT ANALYSIS
OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ego_df <- as.data.frame(ego)
head(ego_df[1:7])
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 
dev.copy(jpeg,filename="Barplot_GEA_Proteins_T1_T0_norm.jpg");
dev.off ();

# #PCA
PCA_raw <- prcomp(t(dat.log), scale. = FALSE) #principal component analysis (PCA)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = meta$condition)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Disease)) +
  ggtitle("PCA plot of the log-transformed proteins") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15,20,30)) + 
  scale_color_manual(values = c("green", "blue", "red"))
dev.copy(jpeg,filename="PCA_proteins.jpg");
dev.off ()
# # t-SNE
tsne(dat.log,
     labels=as.factor(meta$condition),
     controlscale = TRUE,scale=3)
dev.copy(jpeg,filename="t-SNE_proteins_norm2.jpg");
dev.off ()

# GENES
setwd("/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/PBMC_1st_Injection_DIANN_results20072023")
meta <- read_excel("metadata.xlsx")
data<-readr::read_tsv("report.gg_matrix.tsv")
data$Genes%>% duplicated() %>% any()
dat=data[2:55]
rownames(dat)=data$Genes

## ----log2transform1-----------------------------------------------------------
dat.log = log2(dat)
# Create a logical vector indicating which rows have NAs
na_rows = rowSums(is.na(dat.log)) > 0
#remove rows with NAs
dat.log = na.omit(dat.log)
protein_group_filtered = data$Genes[!na_rows]
rownames(dat.log)=protein_group_filtered
## ----boxplot1-----------------------------------------------------------------
boxplot(dat.log,las=2,main="Cellular proteomic data: Genes")
dev.copy(jpeg,filename="Barplot_before_Normalization_Genes.jpg");
dev.off ();
# Here the data is already median centered, we skip the following step. 
dat.log = equalMedianNormalization(dat.log)
boxplot(dat.log,las=2,main="Cellular proteomic data: Genes")
dev.copy(jpeg,filename="Barplot_after_Normalization_Genes.jpg");
dev.off ();
## ----design-------------------------------------------------------------------
x=meta$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)
# The function model.matrix is used to generate the design matrix
design = model.matrix(~0+cond) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
colnames(design) <- c("Stress","T0","T1")
contrast =  makeContrasts(contrasts="T1-Stress",levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)
# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit3, adjust="BH", sort.by ="p", number=Inf, p.value=0.05)
top_genes=na.omit(top_genes)
head(top_genes)
write.csv(top_genes, "DE_Genes_Stress_T1.csv")
addWorksheet(wb, sheetName = 'Stress-T1 (Genes)')
writeData(wb, sheet = 'Stress-T1 (Genes)', x = cbind(rownames(top_genes),top_genes))

EnhancedVolcano(top_genes,
                lab = rownames(top_genes),
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Proteins: Control vs. Heat stroke T1',
                pCutoff = 0.0001,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="Volcano_Genes_Stress_T1.jpg");
dev.off ();

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "SYMBOL")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='SYMBOL'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, BH = "SYMBOL", all = TRUE)
gene_list=na.omit(gene_list)

# GENE ENRICHMENT ANALYSIS
OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_df <- as.data.frame(ego)
head(ego_df[1:7])
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 
dev.copy(jpeg,filename="Barplot_GEA_Genes_Stress_T1.jpg");
dev.off ();

contrast =  makeContrasts(contrasts="T0-Stress",levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)
# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit3, adjust="BH", sort.by ="p", number=Inf, p.value=0.05)
top_genes=na.omit(top_genes)
write.csv(top_genes, "DE_Genes_Stress_T0.csv")
addWorksheet(wb, sheetName = 'Stress-T0 (Genes)')
writeData(wb, sheet = 'Stress-T0 (Genes)', x = cbind(rownames(top_genes),top_genes))
EnhancedVolcano(top_genes,
                lab = rownames(top_genes),
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Proteins: Control vs. Heat stroke T0',
                pCutoff = 0.0001,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="Volcano_Genes_Stress_T0.jpg");
dev.off ();

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "SYMBOL")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='SYMBOL'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, BH = "SYMBOL", all = TRUE)
gene_list=na.omit(gene_list)

# GENE ENRICHMENT ANALYSIS
OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ego_df <- as.data.frame(ego)
head(ego_df[1:7])
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 
dev.copy(jpeg,filename="Barplot_GEA_Genes_Stress_T0.jpg");
dev.off ();

contrast =  makeContrasts(contrasts="T1-T0",levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)
# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit3, adjust="BH", sort.by ="p", number=Inf, p.value=0.05)
top_genes=na.omit(top_genes)
write.csv(top_genes, "DE_Genes_T0_T1.csv")
addWorksheet(wb, sheetName = 'T0-T1 (Genes)')
writeData(wb, sheet = 'T0-T1 (Genes)', x = cbind(rownames(top_genes),top_genes))
# Save the Excel file
saveWorkbook(wb, output_file)
# VOLCANO PLOT
EnhancedVolcano(top_genes,
                lab = rownames(top_genes),
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Proteins: Heat stroke T0 vs. Heat stroke T1',
                pCutoff = 0.0001,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="Volcano_Genes_T0_T1.jpg");
dev.off ();

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "SYMBOL")
top_genes=cbind(top_genes, rownames(top_genes))
colnames(top_genes)[7]='SYMBOL'
head(top_genes)
gene_list <- merge(top_genes, gene_ids, BH = "SYMBOL", all = TRUE)
gene_list=na.omit(gene_list)

# GENE ENRICHMENT ANALYSIS
OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ego_df <- as.data.frame(ego)
head(ego_df[1:7])
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 
dev.copy(jpeg,filename="Barplot_GEA_Genes_T1_T0.jpg");
dev.off ();

# #PCA
PCA_raw <- prcomp(t(dat.log), scale. = FALSE) #principal component analysis (PCA)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = meta$condition)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Disease)) +
  ggtitle("PCA plot of the log-transformed genes") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15,20,30)) + 
  scale_color_manual(values = c("green", "blue", "red"))
dev.copy(jpeg,filename="PCA_genes.jpg");
dev.off ()
# # t-SNE
tsne(dat.log,
     labels=as.factor(meta$condition),
     controlscale = TRUE,scale=3)
dev.copy(jpeg,filename="t-SNE_genes.jpg");
dev.off ()


# PROMOTERS
# setwd("/Users/mariacastillo/Desktop/Proteomics_data/PBMC_1st_Injection_DIANN_results20072023")
# meta <- read_excel("metadata.xlsx")
# data<-readr::read_tsv("report.pr_matrix.tsv")
# data$Protein.Group%>% duplicated() %>% any()
# colnames(data)[1]='UNIPROT'
# data_ = AnnotationDbi::select(org.Hs.eg.db, keys = data$UNIPROT, columns = "ENTREZID", keytype = "UNIPROT")
# data_new <- merge(data_, data, by = "UNIPROT", all = TRUE)
# data_unique <- make_unique(data_new, "ENTREZID", "UNIPROT", delim = ";")
# d_columns <- 7:60
# data_se <- make_se(data_unique, d_columns, meta)
# 
# dat=data[6:59]
# rownames(dat)=data$Genes
# 
# ## ----log2transform1-----------------------------------------------------------
# dat.log = log2(dat)
# # Create a logical vector indicating which rows have NAs
# na_rows = rowSums(is.na(dat.log)) > 0
# #remove rows with NAs
# dat.log = na.omit(dat.log)
# protein_group_filtered = data$Protein.Group[!na_rows]
# rownames(dat.log)=protein_group_filtered
# ## ----boxplot1-----------------------------------------------------------------
# boxplot(dat.log,las=2,main="Cellular proteomic data: Promoters")
# dev.copy(jpeg,filename="Barplot_before_Normalization_Promoters.jpg");
# dev.off ();
# # Here the data is already median centered, we skip the following step. 
# dat.log = equalMedianNormalization(dat.log)
# boxplot(dat.log,las=2,main="Cellular proteomic data: Promoters")
# dev.copy(jpeg,filename="Barplot_after_Normalization_Promoters.jpg");
# dev.off ();
# ## ----design-------------------------------------------------------------------
# x=meta$condition
# x1=str_replace_all(x," ", "_")
# x1=str_replace_all(x1,"-", "_")
# cond = as.factor(x1)
# # The function model.matrix is used to generate the design matrix
# design = model.matrix(~0+cond) # 0 means no intercept for the linear model
# colnames(design) = gsub("cond","",colnames(design))
# colnames(design)[1:3] <- c("Stress","T0","T1")
# contrast =  makeContrasts(contrasts="T1-Stress",levels=design)
# fit1 <- lmFit(dat.log, design)
# fit2 <- contrasts.fit(fit1,contrasts = contrast)
# fit3 <- eBayes(fit2, 0.01)
# DE_genes2 <- decideTests(fit2)
# summary(DE_genes2)
# # Get the top 10 deferentially expressed genes
# top_genes <- topTable(fit3, adjust.method="BY", sort.by ="p", p.value=0.05,number=Inf)
# top_genes=na.omit(top_genes)
# head(top_genes)
# write.csv(top_genes, "DE_Promoters_Stress_T1.csv")
# EnhancedVolcano(top_genes,
#                 lab = rownames(top_genes),
#                 x = 'logFC',
#                 y = 'adj.P.Val',
#                 pCutoff = 0.0001,
#                 FCcutoff = 1,
#                 pointSize = 1.0)
# dev.copy(jpeg,filename="Volcano_Promoters_Stress_T1.jpg");
# dev.off ();
# 
# gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "UNIPROT")
# top_genes=cbind(top_genes, rownames(top_genes))
# colnames(top_genes)[7]='UNIPROT'
# head(top_genes)
# gene_list <- merge(top_genes, gene_ids, BH = "UNIPROT", all = TRUE)
# gene_list=na.omit(gene_list)
# 
# # GENE ENRICHMENT ANALYSIS
# OrgDb='org.Hs.eg.db'
# ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# ego_df <- as.data.frame(ego)
# head(ego_df[1:7])
# #Bar plot of enriched terms.s
# barplot(ego, showCategory=15) 
# dev.copy(jpeg,filename="Barplot_GEA_Promoters_Stress_T1.jpg");
# dev.off ();
# 
# contrast =  makeContrasts(contrasts="T0-Stress",levels=design)
# fit1 <- lmFit(dat.log, design)
# fit2 <- contrasts.fit(fit1,contrasts = contrast)
# fit3 <- eBayes(fit2, 0.01)
# DE_genes2 <- decideTests(fit2)
# summary(DE_genes2)
# # Get the top 10 deferentially expressed genes
# top_genes <- topTable(fit3, adjust="BY", sort.by ="p", p.value=0.05,number=Inf)
# top_genes=na.omit(top_genes)
# write.csv(top_genes, "DE_Promoters_Stress_T0.csv")
# EnhancedVolcano(top_genes,
#                 lab = rownames(top_genes),
#                 x = 'logFC',
#                 y = 'adj.P.Val',
#                 pCutoff = 0.0001,
#                 FCcutoff = 1,
#                 pointSize = 1.0)
# dev.copy(jpeg,filename="Volcano_Promoters_Stress_T0.jpg");
# dev.off ();
# 
# gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "UNIPROT")
# top_genes=cbind(top_genes, rownames(top_genes))
# colnames(top_genes)[7]='UNIPROT'
# head(top_genes)
# gene_list <- merge(top_genes, gene_ids, BH = "UNIPROT", all = TRUE)
# gene_list=na.omit(gene_list)
# 
# # GENE ENRICHMENT ANALYSIS
# OrgDb='org.Hs.eg.db'
# ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# ego_df <- as.data.frame(ego)
# head(ego_df[1:7])
# #Bar plot of enriched terms.s
# barplot(ego, showCategory=15) 
# dev.copy(jpeg,filename="Barplot_GEA_Promoters_Stress_T0.jpg");
# dev.off ();
# 
# contrast =  makeContrasts(contrasts="T1-T0",levels=design)
# fit1 <- lmFit(dat.log, design)
# fit2 <- contrasts.fit(fit1,contrasts = contrast)
# fit3 <- eBayes(fit2, 0.01)
# DE_genes2 <- decideTests(fit2)
# summary(DE_genes2)
# # Get the top 10 deferentially expressed genes
# top_genes <- topTable(fit3, adjust="BY", sort.by="p", p.value=0.05,number=Inf)
# top_genes=na.omit(top_genes)
# write.csv(top_genes, "DE_Promoters_T0_T1.csv")
# # VOLCANO PLOT
# EnhancedVolcano(top_genes,
#                 lab = rownames(top_genes),
#                 x = 'logFC',
#                 y = 'adj.P.Val',
#                 pCutoff = 0.0001,
#                 FCcutoff = 1,
#                 pointSize = 1.0)
# dev.copy(jpeg,filename="Volcano_Promoters_T0_T1.jpg");
# dev.off ();
# 
# gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "ENTREZID", keytype = "UNIPROT")
# top_genes=cbind(top_genes, rownames(top_genes))
# colnames(top_genes)[7]='UNIPROT'
# head(top_genes)
# gene_list <- merge(top_genes, gene_ids, BH = "UNIPROT", all = TRUE)
# gene_list=na.omit(gene_list)
# 
# # GENE ENRICHMENT ANALYSIS
# OrgDb='org.Hs.eg.db'
# ego <- enrichGO(gene_list$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# ego_df <- as.data.frame(ego)
# head(ego_df[1:7])
# #Bar plot of enriched terms.s
# barplot(ego, showCategory=15) 
# dev.copy(jpeg,filename="Barplot_GEA_Promoters_T1_T0.jpg");
# dev.off ();

# #PCA
# 
# vsdata <- vst(dds, nsub=nrow(dds)/2, blind=FALSE)
# plotPCA(vsdata, intgroup="condition") #using the DESEQ2 plotPCA fxn we can
# 
# dev.copy(jpeg,filename="PCA_Promoters.jpg");
# dev.off ()
# 
# # t-SNE
# tsne(dat.log,labels=as.factor(rownames(dat.log)))

library(ggvenn)
library(VennDiagram)
library(openxlsx)
library(topGO)
library(gprofiler2)
library("qpcR")      

transcriptomics_HS_T0 = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_DE_analysis_genes.xlsx', sheet= 'Heat stress - Stroke T0')
proteomic_genes_HS_T0 = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Ressults_proteomics_cell_analysis.xlsx', sheet= 'DE_Genes_Stress_T1')
proteomic_proteins_HS_T0 = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Ressults_proteomics_cell_analysis.xlsx', sheet= 'DE_Proteins_Stress_T1')
Genes_p = AnnotationDbi::select(org.Hs.eg.db, keys = proteomic_proteins_HS_T0$UNIPROT, columns = "SYMBOL", keytype = "UNIPROT")
Genes_pro <- merge(proteomic_proteins_HS_T0, Genes_p, BY = "UNIPROT", all = TRUE)
Genes_proteins=na.omit(Genes_pro$SYMBOL)
methylation_HS_T0 = read.csv('/Users/mariacastillo/Desktop/HEATSTROKE/Methylation/results_final/reports_ALL/differential_methylation_data/diffMethTable_region_cmp1_genes.csv')
methylation_HS_T0=na.omit(methylation_HS_T0)
methylation_HS_T0_F=methylation_HS_T0$symbol

int_trans_meth=list(intersect(transcriptomics_HS_T0$SYMBOL, methylation_HS_T0_F))
int_trans_prote_g=list(intersect(transcriptomics_HS_T0$SYMBOL, proteomic_genes_HS_T0$SYMBOL))
int_trans_prote_p=list(intersect(transcriptomics_HS_T0$SYMBOL, Genes_proteins))
int_prot_meth_g=list(intersect(proteomic_genes_HS_T0$SYMBOL, methylation_HS_T0_F))
int_prot_meth_p=list(intersect(Genes_proteins, methylation_HS_T0_F))

int_trans_meth_trans_prote_p=list(intersect(unlist(int_trans_meth), unlist(int_trans_prote_p)))#transcription,methylation,proteomics(p)
int_trans_meth_trans_prote_g=list(intersect(unlist(int_trans_meth), unlist(int_trans_prote_g)))#transcription,methylation,proteomics(g)

int_trans_meth_prot_meth_p=list(intersect(unlist(int_trans_meth), unlist(int_prot_meth_p)))#transcription,methylation,proteomics(p)
int_trans_meth_prot_meth_g=list(intersect(unlist(int_trans_meth), unlist(int_prot_meth_g)))#transcription,methylation,proteomics(g)

int_trans_prote_prot_meth_p=list(intersect(unlist(int_trans_prote_p), unlist(int_prot_meth_p)))#transcription,methylation,proteomics(p)
int_trans_prote_prot_meth_g=list(intersect(unlist(int_trans_prote_g), unlist(int_prot_meth_g)))#transcription,methylation,proteomics(g)

int_trans_prote_prot_meth_=list(intersect(unlist(int_prot_meth_g), unlist(int_prot_meth_p)))

# Create an Excel file and add sheets with the lists
output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/shared_genes_lists.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'int_trans_meth')
writeData(wb, sheet = 'int_trans_meth', x = int_trans_meth)
addWorksheet(wb, sheetName = 'int_trans_prote(g)')
writeData(wb, sheet = 'int_trans_prote(g)', x = int_trans_prote_g)
addWorksheet(wb, sheetName = 'int_trans_prote(p)')
writeData(wb, sheet = 'int_trans_prote(p)', x = int_trans_prote_p)
addWorksheet(wb, sheetName = 'int_prot(g)_meth')
writeData(wb, sheet = 'int_prot(g)_meth', x = int_prot_meth_g)
addWorksheet(wb, sheetName = 'int_prot(p)_meth')
writeData(wb, sheet = 'int_prot(p)_meth', x = int_prot_meth_p)
addWorksheet(wb, sheetName = 'int_meth_trans_prot(g)')
writeData(wb, sheet = 'int_meth_trans_prot(g)', x = int_trans_prote_prot_meth_g)
addWorksheet(wb, sheetName = 'int_meth_trans_prot(p)')
writeData(wb, sheet = 'int_meth_trans_prot(p)', x = int_trans_prote_prot_meth_p)
addWorksheet(wb, sheetName = 'int_meth_trans_prot(p)_prot(g)')
writeData(wb, sheet = 'int_meth_trans_prot(p)_prot(g)', x = int_trans_prote_prot_meth_)
# Save the Excel file
saveWorkbook(wb, output_file)


GO_int_trans_meth <- gost(query = int_trans_meth, 
                          organism = "hsapiens", ordered_query = TRUE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = FALSE, 
                          user_threshold = 0.05, correction_method = "fdr", 
                          domain_scope = "annotated", custom_bg = NULL, 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_int_trans_prote_g <- gost(query = int_trans_prote_g, 
                             organism = "hsapiens", ordered_query = TRUE, 
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                             measure_underrepresentation = FALSE, evcodes = FALSE, 
                             user_threshold = 0.05, correction_method = "fdr", 
                             domain_scope = "annotated", custom_bg = NULL, 
                             numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_int_trans_prote_p <- gost(query = int_trans_prote_p, 
                             organism = "hsapiens", ordered_query = TRUE, 
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                             measure_underrepresentation = FALSE, evcodes = FALSE, 
                             user_threshold = 0.05, correction_method = "fdr", 
                             domain_scope = "annotated", custom_bg = NULL, 
                             numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_int_prot_meth_g <- gost(query = int_prot_meth_g, 
                           organism = "hsapiens", ordered_query = TRUE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = FALSE, 
                           user_threshold = 0.05, correction_method = "fdr", 
                           domain_scope = "annotated", custom_bg = NULL, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_int_prot_meth_p <- gost(query = int_prot_meth_p, 
                           organism = "hsapiens", ordered_query = TRUE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = FALSE, 
                           user_threshold = 0.05, correction_method = "fdr", 
                           domain_scope = "annotated", custom_bg = NULL, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_int_t_p_p_m_g <- gost(query = int_trans_prote_prot_meth_g, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_int_t_p_p_m_p <- gost(query = int_trans_prote_prot_meth_p, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_int_tran_prot_prot_meth_ <- gost(query = int_trans_prote_prot_meth_, 
                                    organism = "hsapiens", ordered_query = TRUE, 
                                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                                    user_threshold = 0.05, correction_method = "fdr", 
                                    domain_scope = "annotated", custom_bg = NULL, 
                                    numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

# Create an Excel file and add sheets with the lists
output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/shared_genes_GEA.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'GO_int_trans_meth')
writeData(wb, sheet = 'GO_int_trans_meth', x = cbind(GO_int_trans_meth$result$term_id,GO_int_trans_meth$result$term_name))
addWorksheet(wb, sheetName = 'GO_int_trans_prote_g')
writeData(wb, sheet = 'GO_int_trans_prote_g', x = cbind(GO_int_trans_prote_g$result$term_id,GO_int_trans_prote_g$result$term_name))
addWorksheet(wb, sheetName = 'GO_int_trans_prote_p')
writeData(wb, sheet = 'GO_int_trans_prote_p', x = cbind(GO_int_trans_prote_p$result$term_id,GO_int_trans_prote_p$result$term_name))
addWorksheet(wb, sheetName = 'GO_int_prot_meth_g')
writeData(wb, sheet = 'GO_int_prot_meth_g', x = cbind(GO_int_prot_meth_g$result$term_id,GO_int_prot_meth_g$result$term_name))
addWorksheet(wb, sheetName = 'GO_int_prot_meth_p')
writeData(wb, sheet = 'GO_int_prot_meth_p', x = cbind(GO_int_prot_meth_p$result$term_id,GO_int_prot_meth_p$result$term_name))
addWorksheet(wb, sheetName = 'GO_int_t_p_p_m_g')
writeData(wb, sheet = 'GO_int_t_p_p_m_g', x = cbind(GO_int_t_p_p_m_g$result$term_id,GO_int_t_p_p_m_g$result$term_name))
addWorksheet(wb, sheetName = 'GO_int_t_p_p_m_p')
writeData(wb, sheet = 'GO_int_t_p_p_m_p', x = cbind(GO_int_t_p_p_m_p$result$term_id,GO_int_t_p_p_m_p$result$term_name))
addWorksheet(wb, sheetName = 'GO_int_tran_prot_prot_meth_')
writeData(wb, sheet = 'GO_int_tran_prot_prot_meth_', x = cbind(GO_int_tran_prot_prot_meth_$result$term_id,GO_int_tran_prot_prot_meth_$result$term_name))
# Save the Excel file
saveWorkbook(wb, output_file)

transcript_methylation=c(transcriptomics_HS_T0$SYMBOL, methylation_HS_T0_F)
proteomic_g_methylation=c(proteomic_genes_HS_T0$SYMBOL,methylation_HS_T0_F)
transcript_proteomic_g=c(transcriptomics_HS_T0$SYMBOL,proteomic_genes_HS_T0$SYMBOL)
proteomic_p_methylation=c(Genes_proteins, methylation_HS_T0_F)
transcript_proteomic_p=c(transcriptomics_HS_T0$SYMBOL,Genes_proteins)

unique_transcriptomics_g <- setdiff(transcriptomics_HS_T0$SYMBOL,proteomic_g_methylation)
unique_transcriptomics_p <- setdiff(transcriptomics_HS_T0$SYMBOL,proteomic_p_methylation)
unique_methylation_g <- setdiff(methylation_HS_T0_F,transcript_proteomic_g)
unique_methylation_p <- setdiff(methylation_HS_T0_F,transcript_proteomic_p)
unique_proteomics_p <- setdiff(Genes_proteins,transcript_methylation)
unique_proteomics_g <- setdiff(proteomic_genes_HS_T0$SYMBOL,transcript_methylation)


output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/unique_genes_lists.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'unique_transcriptomics_g')
writeData(wb, sheet = 'unique_transcriptomics_g', x = unique_transcriptomics_g)
addWorksheet(wb, sheetName = 'unique_transcriptomics_p')
writeData(wb, sheet = 'unique_transcriptomics_p', x = unique_transcriptomics_p)
addWorksheet(wb, sheetName = 'unique_methylation_g')
writeData(wb, sheet = 'unique_methylation_g', x = unique_methylation_g)
addWorksheet(wb, sheetName = 'unique_methylation_p')
writeData(wb, sheet = 'unique_methylation_p', x = unique_methylation_p)
addWorksheet(wb, sheetName = 'unique_proteomics_p')
writeData(wb, sheet = 'unique_proteomics_p', x = unique_proteomics_p)
addWorksheet(wb, sheetName = 'unique_proteomics_g')
writeData(wb, sheet = 'unique_proteomics_g', x = unique_proteomics_g)
# Save the Excel file
saveWorkbook(wb, output_file)

GO_unique_transcriptomics_g <- gost(query = unique_transcriptomics_g, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_transcriptomics_p <- gost(query = unique_transcriptomics_p, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_methylation_g <- gost(query = unique_methylation_g, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_methylation_p <- gost(query = unique_methylation_p, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_proteomics_p <- gost(query = unique_proteomics_p, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_proteomics_g <- gost(query = unique_proteomics_g, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/unique_genes_GEA.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'unique_transcriptomics_g')
writeData(wb, sheet = 'unique_transcriptomics_g', x = cbind(GO_unique_transcriptomics_g$result$term_id,GO_unique_transcriptomics_g$result$term_name))
addWorksheet(wb, sheetName = 'unique_transcriptomics_p')
writeData(wb, sheet = 'unique_transcriptomics_p', x = cbind(GO_unique_transcriptomics_p$result$term_id,GO_unique_transcriptomics_p$result$term_name))
addWorksheet(wb, sheetName = 'unique_methylation_g')
writeData(wb, sheet = 'unique_methylation_g', x = cbind(GO_unique_methylation_g$result$term_id,GO_unique_methylation_g$result$term_name))
addWorksheet(wb, sheetName = 'unique_methylation_p')
writeData(wb, sheet = 'unique_methylation_p', x = cbind(GO_unique_methylation_p$result$term_id,GO_unique_methylation_p$result$term_name))
addWorksheet(wb, sheetName = 'unique_proteomics_p')
writeData(wb, sheet = 'unique_proteomics_p', x = cbind(GO_unique_proteomics_p$result$term_id,GO_unique_proteomics_p$result$term_name))
addWorksheet(wb, sheetName = 'unique_proteomics_g')
writeData(wb, sheet = 'unique_proteomics_g', x = cbind(GO_unique_proteomics_g$result$term_id,GO_unique_proteomics_g$result$term_name))
# Save the Excel file
saveWorkbook(wb, output_file)

a <- list('Transcriptomis' = transcriptomics_HS_T0$SYMBOL,
          #'Cell proteomics: genes' = proteomic_genes_HS_T0$SYMBOL,
          'Cell proteomics: proteins'= Genes_proteins,
          'Methylation' = methylation_HS_T0_F)

ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Figures/Venn_Diagram_Colors_Trans_Proteo_Meth.jpg");
dev.off ()

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpDisease()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(list(transcriptomics_HS_T0$SYMBOL, proteomic_genes_HS_T0$SYMBOL, methylation_HS_T0_F),category.names = c("Transcriptomis" , "Cell proteomics " , "Methylation"))
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Figures/Venn_Diagram_Trans_Proteo_Meth.jpg");
dev.off ()















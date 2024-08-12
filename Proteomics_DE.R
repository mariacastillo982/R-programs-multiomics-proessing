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
library(SummarizedExperiment)
library(vsn)

# PROTEINS
path_proteomics="/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/PBMC_1st_Injection_DIANN_results20072023/"
setwd(path_proteomics)

meta <- read_excel("metadata.xlsx")
# Read the data from the specified file
data<-readr::read_tsv("report.pg_matrix.tsv")

data$Genes%>% duplicated() %>% any()
# Extract relevant data columns and set rownames
dat <- data[6:ncol(data)]
rownames(dat) <- data$Protein.Group
colnames(dat) <- meta$`Sample ID`

# Log2 transform the data
dat.log <- log2(dat)

# Remove rows with NAs
na_rows <- rowSums(is.na(dat.log)) > 0
dat.log <- na.omit(dat.log)
protein_group_filtered <- data$Protein.Group[!na_rows]
rownames(dat.log) <- protein_group_filtered


meta_Old <- meta[ meta$Age=="Old", ]
meta_Young <- meta[ meta$Age=="Young", ]
meta_Male <- meta[ meta$Gender=="M", ]
meta_Female <- meta[ meta$Gender=="F", ]

rownames(dat)=data$Protein.Group
colnames(dat)=meta$`Sample ID`

gs= factor(meta$condition)
Disease = factor(meta$condition)
gender = factor(meta$Gender)
age = factor(meta$Age)

x=meta$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

process_proteomics_data <- function(dat.log, meta, output_dir) {

  ## ----boxplot1-----------------------------------------------------------------
  # Plot before normalization
  boxplot(dat.log, las=2, main="Cellular proteomic data: Genes")
  dev.copy(jpeg, filename=file.path(output_dir, "Barplot_before_Normalization_Genes.jpg"))
  dev.off()
  
  # Normalize the data (assuming equalMedianNormalization is defined)
  dat.log <- equalMedianNormalization(dat.log)
  
  # Plot after normalization
  boxplot(dat.log, las=2, main="Cellular proteomic data: Normalized")
  dev.copy(jpeg, filename=file.path(output_dir, "Barplot_after_Normalization_Genes.jpg"))
  dev.off()
  
  return(dat.log)
}

dat.log <- process_proteomics_data(dat.log, meta, output_dir)

# Define file paths
output_dir <- "HEATSTROKE/Proteomics_C_data/Results Overall/PROTEINS/"
age_dir <- "HEATSTROKE/Proteomics_C_data/Age_Gender_proteomics/AGE/"
gender_dir <- "HEATSTROKE/Proteomics_C_data/Age_Gender_proteomics/GENDER/"

# Identify the index of the reference protein in the data frame
reference_index <- which(rownames(dat.log) == 'P04406')

# Normalize by the reference protein
normalized_data <- sweep(dat.log, 2, dat.log[reference_index, ], '/')
normalized_data <- as.data.frame(normalized_data)
write.csv(normalized_data, file.path(output_dir, "normalized_data_HKG.csv"))

# Boxplot after normalization
boxplot(normalized_data, las=2, main="Normalization Housekeeping gene")
ggsave(filename=file.path(output_dir, "Barplot_Normalization_Proteins_HKG.jpg"))

# VSN Normalization
dat.log <- equalMedianNormalization(dat.log)
boxplot(dat.log, las=2, main="VSN Normalization: Proteins")
ggsave(filename=file.path(output_dir, "Barplot_after_Normalization_Proteins.jpg"))

# Q-Q plot
value <- as.numeric(unlist(dat.log))
qqnorm(value, pch = 1, frame = FALSE)
qqline(value, col = "steelblue", lwd = 2)
ggsave(filename=file.path(output_dir, "Q-Q_plot_Proteins.jpg"))

# PCA Plotting Function
plot_pca <- function(data, group_col, title_suffix, file_path, x_pc=1, y_pc=2) {
  PCA_raw <- prcomp(t(data), scale. = FALSE)
  percentVar <- round(100 * PCA_raw$sdev^2 / sum(PCA_raw$sdev^2), 1)
  sd_ratio <- sqrt(percentVar[y_pc] / percentVar[x_pc])
  
  dataGG <- data.frame(PC1 = PCA_raw$x[, x_pc], PC2 = PCA_raw$x[, y_pc], group = group_col)
  dataGG_clean <- dataGG %>%
    drop_na() %>%
    group_by(group) %>%
    dplyr::slice(chull(PC1, PC2))
  
  ggplot(dataGG) +
    aes(x = PC1, y = PC2, color = group) +
    geom_point(aes(shape = group)) +
    geom_polygon(data = dataGG_clean,
                 aes(fill = group, color = NULL),
                 alpha = 0.3,
                 show.legend = FALSE) +
    xlab(paste0("PC", x_pc, ", VarExp: ", percentVar[x_pc], "%")) +
    ylab(paste0("PC", y_pc, ", VarExp: ", percentVar[y_pc], "%")) +
    ggtitle(paste("PCA: Proteins -", title_suffix)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = sd_ratio)
  
  #ggsave(filename=file_path)
}

# PCA Plots
plot_pca(dat.log, gs, "Disease (PC1 vs PC2)", file.path(output_dir, "PCA_Proteins_PC1_PC2_Disease.jpg"))
plot_pca(dat.log, age, "Age (PC2 vs PC3)", file.path(age_dir, "PCA_Proteins_PC2_PC3_Age.jpg"))
plot_pca(dat.log, gender, "Gender (PC1 vs PC2)", file.path(gender_dir, "PCA_Proteins_PC1_PC2_Gender.jpg"))

## ----Design-------------------------------------------------------------------
# Clean up the condition names
meta$condition <- gsub("[ -]", "_", meta$condition)

# Create a design matrix
groups <- make.names(c("Stress", "T1", "T0"))
meta$group <- factor(meta$condition, levels = groups)
design <- model.matrix(~ gs + gender + age, data = meta)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design)[1:3] <- c("Stress", "T1", "T0")

# Fit the linear model
fit1 <- lmFit(dat.log, design)

## ----Differential Expression Analysis-----------------------------------------

analyze_DE <- function(contrast_name, fit1, contrast_formula, output_file, plot_title) {
  contrast <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit1, contrast)
  fit3 <- eBayes(fit2, 0.01)
  DE_genes <- decideTests(fit3)
  print(summary(DE_genes))
  
  # Get the top differentially expressed genes
  top_genes <- topTable(fit3, adjust = "BH", sort.by = "p", number = Inf, p.value = 0.05)
  top_genes <- na.omit(top_genes)
  
  # Add gene symbols
  gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "SYMBOL", keytype = "UNIPROT")
  top_genes <- cbind(top_genes, UNIPROT = rownames(top_genes))
  gene_list <- merge(top_genes, gene_ids, by = "UNIPROT", all = TRUE)
  gene_list <- na.omit(gene_list)
  
  # Save to CSV
  write.csv(gene_list, output_file)
  
  # Volcano plot
  EnhancedVolcano(top_genes,
                  lab = rownames(top_genes),
                  x = 'logFC',
                  y = 'adj.P.Val',
                  xlim = c(min(top_genes[['logFC']], na.rm = TRUE) - 0.5, max(top_genes[['logFC']], na.rm = TRUE) + 0.5),
                  ylim = c(0, max(-log10(top_genes[['adj.P.Val']]), na.rm = TRUE) + 5),
                  title = plot_title,
                  pCutoff = 0.0001,
                  FCcutoff = 1.5,
                  pointSize = 1.0)
  dev.copy(jpeg, filename = gsub(".csv", ".jpg", output_file))
  dev.off()
  
  return(gene_list)
}

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(output_dir,"/DE_Proteins_Stress_T0.csv"), "Proteins: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(output_dir,"/DE_Proteins_Stress_T1.csv"), "Proteins: Control vs. Heat stroke T1")
t0_t1 <- analyze_DE("T0-T1", fit1, "T0-T1", paste0(output_dir,"/DE_Proteins_T0_T1.csv"), "Proteins: Heat stroke T0 vs. Heat stroke T1")

## ----Gene Enrichment Analysis-------------------------------------------------
perform_enrichment <- function(gene_list, output_file, plot_title) {
  gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_list$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
  
  OrgDb <- 'org.Hs.eg.db'
  ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  ego_df <- as.data.frame(ego)
  print(head(ego_df[1:7]))
  
  # Bar plot of enriched terms
  barplot(ego, showCategory = 15)
  dev.copy(jpeg, filename = gsub(".csv", "_GEA.jpg", output_file))
  dev.off()
  
  return(ego_df)
}

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(output_dir,"DE_Proteins_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(output_dir,"DE_Proteins_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
ego_t0_t1 <- perform_enrichment(t0_t1, paste0(output_dir,"DE_Proteins_T0_T1.csv"), "Gene Enrichment: T0 vs. T1")

# t-SNE
tsne(dat.log,
     labels=as.factor(meta$condition),
     controlscale = TRUE,scale=3)
dev.copy(jpeg,filename=paste0(output_dir,"/t-SNE_proteins_norm2.jpg"));
dev.off ();

## ----AGE-GENDER---------------------------------------------------------------

dat_Old = dat.log[,as.character(meta_Old$`Sample ID`)]
dat_Young = dat.log[,as.character(meta_Young$`Sample ID`)]

dat_Male = dat.log[,as.character(meta_Male$`Sample ID`)]
dat_Female = dat.log[,as.character(meta_Female$`Sample ID`)]


# ----AGE----------------------
#______________________________YOUNG_______________________________
gs= factor(meta_Young$condition)
gender = factor(meta_Young$Gender)
x=meta_Young$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))
# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+gender) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T0","T1")
# Fit the linear model
fit1 <- lmFit(dat_Young, design)

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(age_dir,"YOUNG/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(age_dir,"YOUNG/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")
#t0_t1 <- analyze_DE("T0-T1", fit1, "T0-T1", paste0(age_dir,"YOUNG/DE_Genes_T0_T1.csv"), "Genes: Heat stroke T0 vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(age_dir,"YOUNG/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(age_dir,"YOUNG/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
                                                                                                                                            ego_t0_t1 <- perform_enrichment(t0_t1, paste0(age_dir,"YOUNG/DE_Genes_T0_T1.csv", "Gene Enrichment: T0 vs. T1")
#______________________________OLD_______________________________
gs= factor(meta_Old$condition)
gender = factor(meta_Old$Gender)
x=meta_Old$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+gender) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T1","T0")
# Fit the linear model
fit1 <- lmFit(dat_Old, design)

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(age_dir,"OLD/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(age_dir,"OLD/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(age_dir,"OLD/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(age_dir,"OLD/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")

                                                                                                                                                          
# ----GENDER----------------------

#______________________________MALE_______________________________
gs= factor(meta_Male$condition)
age = factor(meta_Male$Age)

x=meta_Male$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+age) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T1","T0")
# Fit the linear model
fit1 <- lmFit(dat_Male, design)
# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(gender_dir,"MALE/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(gender_dir,"MALE/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(gender_dir,"MALE/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(gender_dir,"MALE/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")

#______________________________FEMALE_______________________________
gs= factor(meta_Female$condition)
age = factor(meta_Female$Age)

x=meta_Female$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+age) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T1","T0")
# Fit the linear model
fit1 <- lmFit(dat_Female, design)
# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(gender_dir,"FEMALE/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(gender_dir,"FEMALE/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(gender_dir,"FEMALE/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(gender_dir,"FEMALE/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
    
# ______________________________________GENES___________________________________

readr::read_tsv("report.gg_matrix.tsv")
data$Genes%>% duplicated() %>% any()
dat=data[2:55]
rownames(dat)=data$Genes
colnames(dat)=meta$Sample ID
output_dir="/path/HEATSTROKE/Proteomics_C_data/Results Overall/GENES"

dat.log <- process_proteomics_data(dat, meta, output_dir)

## ----design-------------------------------------------------------------------
# The function model.matrix is used to generate the design matrix
gs= factor(meta$condition)
Disease = factor(meta$condition)
gender = factor(meta$Gender)
age = factor(meta$Age)
x=meta$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

design = model.matrix(~cond+gender+age) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T0","T1")
rownames(design)=cond
# Fit the linear model
fit1 <- lmFit(dat.log, design)

# PCA Plots
plot_pca(dat.log, gs, "Disease (PC1 vs PC2)", file.path(output_dir, "PCA_Genes_PC1_PC2_Disease.jpg"))
plot_pca(dat.log, age, "Age (PC2 vs PC3)", file.path(age_dir, "PCA_Genes_PC2_PC3_Age.jpg"))
plot_pca(dat.log, gender, "Gender (PC1 vs PC2)", file.path(gender_dir, "PCA_Genes_PC1_PC2_Gender.jpg"))

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(output_dir,"GENES/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(output_dir,"GENES/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(output_dir,"GENES/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(output_dir,"GENES/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")

# # t-SNE
tsne(dat.log,
     labels=as.factor(meta$condition),
     controlscale = TRUE,scale=3)
dev.copy(jpeg,filename=paste0(output_dir,"/HEATSTROKE/Proteomics_C_data/Results Overall/GENES/t-SNE_genes.jpg"));
dev.off ()

a <- list('DE T0' = rownames(top_genes_t0),
          'DE T1' = rownames(top_genes_t1))
venn <- ggvenn(a)
# Add a title to the Venn diagram
venn <- venn + ggtitle("DE Proteins")
# Print the Venn diagram
print(venn)
dev.copy(jpeg,filename="path/HEATSTROKE/Proteomics_C_data/Results Overall/GENES/DE_genes_T0_T1.jpg");
dev.off ();

# PROMOTERS

data<-readr::read_tsv("report.pr_matrix.tsv")
data$Protein.Group%>% duplicated() %>% any()
colnames(data)[1]='UNIPROT'
data_ = AnnotationDbi::select(org.Hs.eg.db, keys = data$UNIPROT, columns = "ENTREZID", keytype = "UNIPROT")
data_new <- merge(data_, data, by = "UNIPROT", all = TRUE)
data_unique <- make_unique(data_new, "ENTREZID", "UNIPROT", delim = ";")
d_columns <- 7:60
data_se <- make_se(data_unique, d_columns, meta)

dat=data[6:59]
rownames(dat)=data$Genes

dat.log <- process_proteomics_data(dat, meta, output_dir)

# PCA Plots
plot_pca(dat.log, gs, "Disease (PC1 vs PC2)", file.path(output_dir, "PCA_Promoters_PC1_PC2_Disease.jpg"))
plot_pca(dat.log, age, "Age (PC2 vs PC3)", file.path(age_dir, "PCA_Promoters_PC2_PC3_Age.jpg"))
plot_pca(dat.log, gender, "Gender (PC1 vs PC2)", file.path(gender_dir, "PCA_Promoters_PC1_PC2_Gender.jpg"))

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(output_dir,"Promoters/DE_Promoters_Stress_T0.csv"), "Promoters: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(output_dir,"Promoters/DE_Promoters_Stress_T1.csv"), "Promoters: Control vs. Heat stroke T1")
t0_t1 <- analyze_DE("T0-T1", fit1, "T0-T1", paste0(output_dir,"Promoters/DE_Promoters_T0_T1.csv"), "Promoters: Heat stroke T0 vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(output_dir,"Promoters/DE_Promoters_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(output_dir,"Promoters/DE_Promoters_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
ego_t0_t1 <- perform_enrichment(t0_t1, paste0(output_dir,"Promoters/DE_Promoters_T0_T1.csv"), "Gene Enrichment: T0 vs. T1")

# t-SNE
tsne(dat.log,labels=as.factor(rownames(dat.log)))## ----LoadpackDisease--------------------------------------------------------------
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
library(SummarizedExperiment)
library(vsn)

# PROTEINS
path_proteomics="/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/PBMC_1st_Injection_DIANN_results20072023/"
setwd(path_proteomics)

meta <- read_excel("metadata.xlsx")
# Read the data from the specified file
data<-readr::read_tsv("report.pg_matrix.tsv")

data$Genes%>% duplicated() %>% any()
# Extract relevant data columns and set rownames
dat <- data[6:ncol(data)]
rownames(dat) <- data$Protein.Group
colnames(dat) <- meta$`Sample ID`

# Log2 transform the data
dat.log <- log2(dat)

# Remove rows with NAs
na_rows <- rowSums(is.na(dat.log)) > 0
dat.log <- na.omit(dat.log)
protein_group_filtered <- data$Protein.Group[!na_rows]
rownames(dat.log) <- protein_group_filtered


meta_Old <- meta[ meta$Age=="Old", ]
meta_Young <- meta[ meta$Age=="Young", ]
meta_Male <- meta[ meta$Gender=="M", ]
meta_Female <- meta[ meta$Gender=="F", ]
#dat_Old = dat[,meta_Old$label]
#rownames(dat_Old)=data$Protein.Group
#colnames(dat_Old)=meta_Old$`Sample ID`
#dat_Young = dat[,meta_Young$label]
#rownames(dat_Young)=data$Protein.Group
#colnames(dat_Young)=meta_Young$`Sample ID`

rownames(dat)=data$Protein.Group
colnames(dat)=meta$`Sample ID`

gs= factor(meta$condition)
Disease = factor(meta$condition)
gender = factor(meta$Gender)
age = factor(meta$Age)

x=meta$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

process_proteomics_data <- function(dat.log, meta, output_dir) {

  ## ----boxplot1-----------------------------------------------------------------
  # Plot before normalization
  boxplot(dat.log, las=2, main="Cellular proteomic data: Genes")
  dev.copy(jpeg, filename=file.path(output_dir, "Barplot_before_Normalization_Genes.jpg"))
  dev.off()
  
  # Normalize the data (assuming equalMedianNormalization is defined)
  dat.log <- equalMedianNormalization(dat.log)
  
  # Plot after normalization
  boxplot(dat.log, las=2, main="Cellular proteomic data: Normalized")
  dev.copy(jpeg, filename=file.path(output_dir, "Barplot_after_Normalization_Genes.jpg"))
  dev.off()
  
  return(dat.log)
}

output_dir="/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Results Overall/PROTEINS/"
dat.log <- process_proteomics_data(dat.log, meta, output_dir)

# Define file paths
output_dir <- "/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Results Overall/PROTEINS/"
age_dir <- "/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Age_Gender_proteomics/AGE/"
gender_dir <- "/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Age_Gender_proteomics/GENDER/"

# Identify the index of the reference protein in the data frame
reference_index <- which(rownames(dat.log) == 'P04406')

# Normalize by the reference protein
normalized_data <- sweep(dat.log, 2, dat.log[reference_index, ], '/')
normalized_data <- as.data.frame(normalized_data)
write.csv(normalized_data, file.path(output_dir, "normalized_data_HKG.csv"))

# Boxplot after normalization
boxplot(normalized_data, las=2, main="Normalization Housekeeping gene")
ggsave(filename=file.path(output_dir, "Barplot_Normalization_Proteins_HKG.jpg"))

# VSN Normalization
dat.log <- equalMedianNormalization(dat.log)
boxplot(dat.log, las=2, main="VSN Normalization: Proteins")
ggsave(filename=file.path(output_dir, "Barplot_after_Normalization_Proteins.jpg"))

# Q-Q plot
value <- as.numeric(unlist(dat.log))
qqnorm(value, pch = 1, frame = FALSE)
qqline(value, col = "steelblue", lwd = 2)
ggsave(filename=file.path(output_dir, "Q-Q_plot_Proteins.jpg"))

# PCA Plotting Function
plot_pca <- function(data, group_col, title_suffix, file_path, x_pc=1, y_pc=2) {
  PCA_raw <- prcomp(t(data), scale. = FALSE)
  percentVar <- round(100 * PCA_raw$sdev^2 / sum(PCA_raw$sdev^2), 1)
  sd_ratio <- sqrt(percentVar[y_pc] / percentVar[x_pc])
  
  dataGG <- data.frame(PC1 = PCA_raw$x[, x_pc], PC2 = PCA_raw$x[, y_pc], group = group_col)
  dataGG_clean <- dataGG %>%
    drop_na() %>%
    group_by(group) %>%
    dplyr::slice(chull(PC1, PC2))
  
  ggplot(dataGG) +
    aes(x = PC1, y = PC2, color = group) +
    geom_point(aes(shape = group)) +
    geom_polygon(data = dataGG_clean,
                 aes(fill = group, color = NULL),
                 alpha = 0.3,
                 show.legend = FALSE) +
    xlab(paste0("PC", x_pc, ", VarExp: ", percentVar[x_pc], "%")) +
    ylab(paste0("PC", y_pc, ", VarExp: ", percentVar[y_pc], "%")) +
    ggtitle(paste("PCA: Proteins -", title_suffix)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = sd_ratio)
  
  #ggsave(filename=file_path)
}

# PCA Plots
plot_pca(dat.log, gs, "Disease (PC1 vs PC2)", file.path(output_dir, "PCA_Proteins_PC1_PC2_Disease.jpg"))
plot_pca(dat.log, age, "Age (PC2 vs PC3)", file.path(age_dir, "PCA_Proteins_PC2_PC3_Age.jpg"))
plot_pca(dat.log, gender, "Gender (PC1 vs PC2)", file.path(gender_dir, "PCA_Proteins_PC1_PC2_Gender.jpg"))

## ----Design-------------------------------------------------------------------
# Clean up the condition names
meta$condition <- gsub("[ -]", "_", meta$condition)

# Create a design matrix
groups <- make.names(c("Stress", "T1", "T0"))
meta$group <- factor(meta$condition, levels = groups)
design <- model.matrix(~ gs + gender + age, data = meta)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design)[1:3] <- c("Stress", "T1", "T0")

# Fit the linear model
fit1 <- lmFit(dat.log, design)

## ----Differential Expression Analysis-----------------------------------------

analyze_DE <- function(contrast_name, fit1, contrast_formula, output_file, plot_title) {
  contrast <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit1, contrast)
  fit3 <- eBayes(fit2, 0.01)
  DE_genes <- decideTests(fit3)
  print(summary(DE_genes))
  
  # Get the top differentially expressed genes
  top_genes <- topTable(fit3, adjust = "BH", sort.by = "p", number = Inf, p.value = 0.05)
  top_genes <- na.omit(top_genes)
  
  # Add gene symbols
  gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(top_genes), columns = "SYMBOL", keytype = "UNIPROT")
  top_genes <- cbind(top_genes, UNIPROT = rownames(top_genes))
  gene_list <- merge(top_genes, gene_ids, by = "UNIPROT", all = TRUE)
  gene_list <- na.omit(gene_list)
  
  # Save to CSV
  write.csv(gene_list, output_file)
  
  # Volcano plot
  EnhancedVolcano(top_genes,
                  lab = rownames(top_genes),
                  x = 'logFC',
                  y = 'adj.P.Val',
                  xlim = c(min(top_genes[['logFC']], na.rm = TRUE) - 0.5, max(top_genes[['logFC']], na.rm = TRUE) + 0.5),
                  ylim = c(0, max(-log10(top_genes[['adj.P.Val']]), na.rm = TRUE) + 5),
                  title = plot_title,
                  pCutoff = 0.0001,
                  FCcutoff = 1.5,
                  pointSize = 1.0)
  dev.copy(jpeg, filename = gsub(".csv", ".jpg", output_file))
  dev.off()
  
  return(gene_list)
}

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(output_dir,"/DE_Proteins_Stress_T0.csv"), "Proteins: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(output_dir,"/DE_Proteins_Stress_T1.csv"), "Proteins: Control vs. Heat stroke T1")
t0_t1 <- analyze_DE("T0-T1", fit1, "T0-T1", paste0(output_dir,"/DE_Proteins_T0_T1.csv"), "Proteins: Heat stroke T0 vs. Heat stroke T1")

## ----Gene Enrichment Analysis-------------------------------------------------
perform_enrichment <- function(gene_list, output_file, plot_title) {
  gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_list$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
  
  OrgDb <- 'org.Hs.eg.db'
  ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  ego_df <- as.data.frame(ego)
  print(head(ego_df[1:7]))
  
  # Bar plot of enriched terms
  barplot(ego, showCategory = 15)
  dev.copy(jpeg, filename = gsub(".csv", "_GEA.jpg", output_file))
  dev.off()
  
  return(ego_df)
}

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(output_dir,"DE_Proteins_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(output_dir,"DE_Proteins_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
ego_t0_t1 <- perform_enrichment(t0_t1, paste0(output_dir,"DE_Proteins_T0_T1.csv"), "Gene Enrichment: T0 vs. T1")

# t-SNE
tsne(dat.log,
     labels=as.factor(meta$condition),
     controlscale = TRUE,scale=3)
dev.copy(jpeg,filename=paste0(output_dir,"/t-SNE_proteins_norm2.jpg"));
dev.off ();

## ----AGE-GENDER---------------------------------------------------------------

dat_Old = dat.log[,as.character(meta_Old$`Sample ID`)]
dat_Young = dat.log[,as.character(meta_Young$`Sample ID`)]

dat_Male = dat.log[,as.character(meta_Male$`Sample ID`)]
dat_Female = dat.log[,as.character(meta_Female$`Sample ID`)]


# ----AGE----------------------
#______________________________YOUNG_______________________________
gs= factor(meta_Young$condition)
gender = factor(meta_Young$Gender)
x=meta_Young$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))
# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+gender) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T0","T1")
# Fit the linear model
fit1 <- lmFit(dat_Young, design)

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(age_dir,"YOUNG/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(age_dir,"YOUNG/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")
#t0_t1 <- analyze_DE("T0-T1", fit1, "T0-T1", paste0(age_dir,"YOUNG/DE_Genes_T0_T1.csv"), "Genes: Heat stroke T0 vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(age_dir,"YOUNG/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(age_dir,"YOUNG/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
                                                                                                                                            ego_t0_t1 <- perform_enrichment(t0_t1, paste0(age_dir,"YOUNG/DE_Genes_T0_T1.csv", "Gene Enrichment: T0 vs. T1")
#______________________________OLD_______________________________
gs= factor(meta_Old$condition)
gender = factor(meta_Old$Gender)
x=meta_Old$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+gender) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T1","T0")
# Fit the linear model
fit1 <- lmFit(dat_Old, design)

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(age_dir,"OLD/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(age_dir,"OLD/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(age_dir,"OLD/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(age_dir,"OLD/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")

                                                                                                                                                          
# ----GENDER----------------------

#______________________________MALE_______________________________
gs= factor(meta_Male$condition)
age = factor(meta_Male$Age)

x=meta_Male$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+age) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T1","T0")
# Fit the linear model
fit1 <- lmFit(dat_Male, design)
# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(gender_dir,"MALE/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(gender_dir,"MALE/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(gender_dir,"MALE/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(gender_dir,"MALE/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")

#______________________________FEMALE_______________________________
gs= factor(meta_Female$condition)
age = factor(meta_Female$Age)

x=meta_Female$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~cond+age) # 0 means no intercept for the linear model
rownames(design) = cond
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T1","T0")
# Fit the linear model
fit1 <- lmFit(dat_Female, design)
# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(gender_dir,"FEMALE/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(gender_dir,"FEMALE/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(gender_dir,"FEMALE/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(gender_dir,"FEMALE/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
    
# ______________________________________GENES___________________________________

readr::read_tsv("report.gg_matrix.tsv")
data$Genes%>% duplicated() %>% any()
dat=data[2:55]
rownames(dat)=data$Genes
colnames(dat)=meta$Sample ID
output_dir="/path/HEATSTROKE/Proteomics_C_data/Results Overall/GENES"

dat.log <- process_proteomics_data(dat, meta, output_dir)

## ----design-------------------------------------------------------------------
# The function model.matrix is used to generate the design matrix
gs= factor(meta$condition)
Disease = factor(meta$condition)
gender = factor(meta$Gender)
age = factor(meta$Age)
x=meta$condition
x1=str_replace_all(x," ", "_")
x1=str_replace_all(x1,"-", "_")
cond = as.factor(x1)

design = model.matrix(~cond+gender+age) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
colnames(design)[1:3] <- c("Stress","T0","T1")
rownames(design)=cond
# Fit the linear model
fit1 <- lmFit(dat.log, design)

# PCA Plots
plot_pca(dat.log, gs, "Disease (PC1 vs PC2)", file.path(output_dir, "PCA_Genes_PC1_PC2_Disease.jpg"))
plot_pca(dat.log, age, "Age (PC2 vs PC3)", file.path(age_dir, "PCA_Genes_PC2_PC3_Age.jpg"))
plot_pca(dat.log, gender, "Gender (PC1 vs PC2)", file.path(gender_dir, "PCA_Genes_PC1_PC2_Gender.jpg"))

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(output_dir,"GENES/DE_Genes_Stress_T0.csv"), "Genes: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(output_dir,"GENES/DE_Genes_Stress_T1.csv"), "Genes: Control vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(output_dir,"GENES/DE_Genes_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(output_dir,"GENES/DE_Genes_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")

# # t-SNE
tsne(dat.log,
     labels=as.factor(meta$condition),
     controlscale = TRUE,scale=3)
dev.copy(jpeg,filename=paste0(output_dir,"/HEATSTROKE/Proteomics_C_data/Results Overall/GENES/t-SNE_genes.jpg"));
dev.off ()

a <- list('DE T0' = rownames(top_genes_t0),
          'DE T1' = rownames(top_genes_t1))
venn <- ggvenn(a)
# Add a title to the Venn diagram
venn <- venn + ggtitle("DE Proteins")
# Print the Venn diagram
print(venn)
dev.copy(jpeg,filename="path/HEATSTROKE/Proteomics_C_data/Results Overall/GENES/DE_genes_T0_T1.jpg");
dev.off ();

# PROMOTERS

data<-readr::read_tsv("report.pr_matrix.tsv")
data$Protein.Group%>% duplicated() %>% any()
colnames(data)[1]='UNIPROT'
data_ = AnnotationDbi::select(org.Hs.eg.db, keys = data$UNIPROT, columns = "ENTREZID", keytype = "UNIPROT")
data_new <- merge(data_, data, by = "UNIPROT", all = TRUE)
data_unique <- make_unique(data_new, "ENTREZID", "UNIPROT", delim = ";")
d_columns <- 7:60
data_se <- make_se(data_unique, d_columns, meta)

dat=data[6:59]
rownames(dat)=data$Genes

dat.log <- process_proteomics_data(dat, meta, output_dir)

# PCA Plots
plot_pca(dat.log, gs, "Disease (PC1 vs PC2)", file.path(output_dir, "PCA_Promoters_PC1_PC2_Disease.jpg"))
plot_pca(dat.log, age, "Age (PC2 vs PC3)", file.path(age_dir, "PCA_Promoters_PC2_PC3_Age.jpg"))
plot_pca(dat.log, gender, "Gender (PC1 vs PC2)", file.path(gender_dir, "PCA_Promoters_PC1_PC2_Gender.jpg"))

# Run the analysis for each comparison
stress_t0 <- analyze_DE("Stress-T0", fit1, "Stress-T0", paste0(output_dir,"Promoters/DE_Promoters_Stress_T0.csv"), "Promoters: Control vs. Heat stroke T0")
stress_t1 <- analyze_DE("Stress-T1", fit1, "Stress-T1", paste0(output_dir,"Promoters/DE_Promoters_Stress_T1.csv"), "Promoters: Control vs. Heat stroke T1")
t0_t1 <- analyze_DE("T0-T1", fit1, "T0-T1", paste0(output_dir,"Promoters/DE_Promoters_T0_T1.csv"), "Promoters: Heat stroke T0 vs. Heat stroke T1")

# Perform enrichment analysis for each comparison
ego_stress_t0 <- perform_enrichment(stress_t0, paste0(output_dir,"Promoters/DE_Promoters_Stress_T0.csv"), "Gene Enrichment: Stress vs. T0")
ego_stress_t1 <- perform_enrichment(stress_t1, paste0(output_dir,"Promoters/DE_Promoters_Stress_T1.csv"), "Gene Enrichment: Stress vs. T1")
ego_t0_t1 <- perform_enrichment(t0_t1, paste0(output_dir,"Promoters/DE_Promoters_T0_T1.csv"), "Gene Enrichment: T0 vs. T1")

# t-SNE
tsne(dat.log,labels=as.factor(rownames(dat.log)))


#Load libraries 
library(affy)
library(affycoretools)
library("readxl")
library(mixtools)
library(clariomdhumantranscriptcluster.db)
library(dplyr)
library(limma)
library(topGO)
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)
library("GEOquery")
library(genefilter)
library(clusterProfiler)
library(clariomdhumanhsrefseq.db)
library(clariomdhumanhsrefseqcdf)
library(clariomdhumanhsrefseqprobe)
library(pd.clariomdhuman.hs.refseq)
library(EnhancedVolcano)
library(ggvenn)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggforce)
library("bladderbatch")
library(nortest)

#Set path where the files are located
setwd("/Users/mariacastillo/Desktop/HEATSTROKE/DATA CEL FILES/ALL")

SDRF <- read_excel("metadata.xlsx")
celFiles <- SDRF$File.name

newname<-list("Sample.number","File.name",
              "Group.number",
              "Group.name", 
              "Time.point", "Gender", "Age","Age.category","Batch")

rownames(SDRF) <- SDRF$File.name
SDRF <- AnnotatedDataFrame(SDRF)

# Reading in .cel files
raw_data<-ReadAffy(verbose = FALSE, filenames =celFiles, phenoData = SDRF, cdfname="clariomdhumanhsrefseq")
stopifnot(validObject(raw_data))
head(Biobase::pData(raw_data))
Biobase::pData(raw_data) <- Biobase::pData(raw_data)
head(Biobase::pData(raw_data))

eset<-affy::rma(raw_data)

#RLE
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(eset)))

RLE_data <- sweep(Biobase::exprs(eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))
dev.copy(jpeg,filename="./PCA/RLE_batch.jpg");
dev.off ();

#HEATMAP
#eset_ <- eset[ ,eset@phenoData@data[["Time.point"]]!="Control1"]
exp_palmieri <- Biobase::exprs(eset)

Gender=factor(eset_@phenoData@data[["Gender"]])
gs <- factor(eset_@phenoData@data[["Time.point"]])
Age=factor(eset_@phenoData@data[["Age.category"]])
batch=factor(eset_@phenoData@data[["Batch"]])

annotation_for_heatmap <- data.frame(Disease = gs, Gender = Gender, Age = Age)

row.names(annotation_for_heatmap) <- eset_@phenoData@data[["Sample number"]]#row.names(pData(eset))
dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"))

rownames(dists) <- eset_@phenoData@data$`Sample number`#eset_@phenoData@data[["Sample number"]]#row.names(pData(eset_))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

pheatmap(dists, col = (hmcol), 
         annotation_row = annotation_for_heatmap,
         legend = TRUE, 
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering Heatmap")

dev.copy(jpeg,filename="./PCA/Heatmap.jpg");
dev.off ();

#PCA
exprs=exprs(eset)
PCA_raw <- prcomp(t(exprs), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

# Calculate the convex hull for each group
dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./PCA/PCA.jpg")
dev.off()

ad.test(exprs)
#create Q-Q plot for both datasets
qqnorm(exprs)
qqline(exprs)
dev.copy(jpeg,filename="./PCA/Q-Q plot.jpg");
dev.off ();

#Filter gene set
variance <- apply(exprs, 1, var)
variance <- na.omit(variance)
model <- normalmixEM(variance, k = 2,  fast=TRUE)
means=model$mu
pos_prob=model$posterior
if (means[1] > means[2]) {
  selectedProbs <- pos_prob[, 1]  # Select column 1
} else {
  selectedProbs <- pos_prob[, 2]  # Select column 2
}
indexes <- which(selectedProbs > 0.5)
palmieri_manfiltered <- eset[indexes, ]

# Annotation of the transcript clusters
anno_palmieri <- AnnotationDbi::select(clariomdhumanhsrefseq.db,
                                       keys = (featureNames(palmieri_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))

# Removing multiple mappings
anno_grouped <- group_by(anno_palmieri, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)
probe_stats <- anno_filtered 
nrow(probe_stats)
ids_to_exclude <- (featureNames(palmieri_manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exclude)
palmieri_final <- subset(palmieri_manfiltered, !ids_to_exclude)
validObject(palmieri_final)
head(anno_palmieri)

fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))
fData(palmieri_final) <- left_join(fData(palmieri_final),anno_palmieri)
#rownames(fData(palmieri_final)) <- anno_palmieri$PROBEID
rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID 
validObject(palmieri_final)

exprs <- exprs(palmieri_final)
rownames(exprs)=c(palmieri_final@featureData@data[["SYMBOL"]])
colnames(exprs)=palmieri_final@phenoData@data[["Time.point"]]
write.csv(exprs, file = "filtered_normalized_exprs.csv",row.names = TRUE, col.names = TRUE)

#remove outlier
#palmieri_final=palmieri_final[,palmieri_final@phenoData@data$File.name!="200-18-000110.CEL"]

# Function to aggregate strings or calculate mean for non-characters
aggregate_strings <- function(x) {
  if (is.character(x)) {
    return(ifelse(length(x) > 1, x[1], x))
  } else {
    return(mean(x))
  }
}

# Filter the palmieri_final_ data based on specific time points
filtered_data <- palmieri_final[, palmieri_final@phenoData@data$"Time.point" %in% c("Control", "Heat Stroke T1", "Heat Stroke T0")]
Subject <- factor(filtered_data@phenoData@data$Group.number)
gs <- factor(filtered_data@phenoData@data$Time.point)
gender <- factor(filtered_data@phenoData@data$Gender)
age <- factor(filtered_data@phenoData@data$Age.category)
batch <- factor(filtered_data@phenoData@data$Batch)

# Assign levels to gs and update palmieri_final_ with group
levels(gs) <- make.names(c("Stress", "T1", "T0"))
filtered_data$group <- gs

# Design matrix for linear model
design_matrix <- model.matrix(~0 + gs + age + gender)
rownames(design_matrix) <- filtered_data@phenoData@data$Time.point
colnames(design_matrix)[1:3] <- c("Stress", "T0", "T1")

# Remove batch effects
norm_cpn <- removeBatchEffect(filtered_data, batch = batch, design = design_matrix)
rownames(norm_cpn) <- filtered_data@featureData@data$SYMBOL

# Function to analyze differential expression for a specific contrast
analyze_contrast <- function(design, contrast, filename_prefix) {
  # Fit linear model and compute statistics
  fit <- lmFit(norm_cpn, design)
  cont.matrix <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  
  # Top table and aggregation
  top_genes <- topTable(fit2, adjust = "BH", sort.by = "B", p.value = 0.05, number = Inf)
  top_genes <- na.omit(top_genes)
  colnames(top_genes)[1] <- "SYMBOL"
  top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
  write.csv(top_genes_mean_table, paste0(filename_prefix, "_top_genes.csv"))
  
  # Add ENTREZ IDs
  gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
  gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)
  write.csv(gene_list, paste0(filename_prefix, "_gene_list.csv"))
  
  # Compare with reference genes
  Ref <- read_excel("Heat_stroke_genes.xlsx")
  merged_table <- merge(Ref, top_genes_mean_table, by = "SYMBOL")
  write.csv(merged_table, paste0(filename_prefix, "_shared_genes.csv"))
  
  # Plot volcano plot
  top_genes_mean_table$logFC <- as.numeric(unlist(top_genes_mean_table$logFC))
  top_genes_mean_table$adj.P.Val <- as.numeric(unlist(top_genes_mean_table$adj.P.Val))
  EnhancedVolcano(top_genes_mean_table,
                  lab = top_genes_mean_table$SYMBOL,
                  x = 'logFC',
                  y = 'adj.P.Val',
                  title = paste('Control vs', contrast),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 1.0)
  dev.copy(jpeg, filename = paste0("Volcano_", filename_prefix, ".jpg"))
  dev.off()
  
  # Enrichment analysis
  OrgDb <- 'org.Hs.eg.db'
  ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  barplot(ego, showCategory = 20)
  dev.copy(jpeg, filename = paste0("Barplot_", filename_prefix, ".jpg"))
  dev.off()
  
  return(top_genes_mean_table)
}

# Analyze for T0
top_genes_T0 <- analyze_contrast(design_matrix, "T0", "Stress_T0_adj_batch_limma")

# Analyze for T1
top_genes_T1 <- analyze_contrast(design_matrix, "T1", "Stress_T1_adj_batch_limma")

# Venn diagram of shared genes
venn_data <- list('DE T0' = top_genes_T0$SYMBOL, 'DE T1' = top_genes_T1$SYMBOL)
venn <- ggvenn(venn_data) + ggtitle("Shared DE Genes between T0 and T1")
print(venn)
dev.copy(jpeg, filename = "DE_genes_T0_T1.jpg");
dev.off();

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
dev.copy(jpeg,filename="RLE_batch.jpg");
dev.off ();

#Heatmap
eset_ <- eset[ ,eset@phenoData@data[["Time.point"]]!="Control1"]
exp_palmieri <- Biobase::exprs(eset_)

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
         legend_labels = (c("small distance", "large distance")))
         #main = "Clustering heatmap for the calibrated samples")

dev.copy(jpeg,filename="Heatmap_gender.jpg");
dev.off ();

PCA_raw <- prcomp(t(y2), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)
library(ggforce)
dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG[dataGG$Disease=='Control',]) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "PCA_batch_PC2_PC3_control_batch_adj.jpg")
dev.off()
  
ggplot(dataGG[dataGG=='Control',]) +
  aes(x = PC1, y = PC2, color = Disease) +  # Color by Disease
  geom_polygon(data = dataGG_clean,
               aes(),  # No need for fill or color
               alpha = 0.3,
               show.legend = TRUE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)  # Assuming sd_ratio is defined
  
dev.copy(jpeg, filename = "PCA_batch_PC1_PC2_control_batch_adj.jpg")
dev.off()

exp_palmieri=exprs(palmieri_final)
y2 <- limma::removeBatchEffect(exp_palmieri, batch = batch)
boxplot(as.data.frame(exp_palmieri),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")

BiocManager::install("bladderbatch")
library("bladderbatch")

pheno = pData(eset)
edata = exprs(eset)
batch = pheno$Batch
gs=factor(pheno$Batch)
age=factor(pheno$Age.category)
mod = model.matrix(~0+gs+age, data=pheno)
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

ggplot(dataGG)+aes(x=PC1,y=PC2,color=Disease)+
  geom_point(aes(shape=Age))+
  stat_ellipse()+stat_ellipse(level=0.8)+
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

dev.copy(jpeg, filename = "PCA_option2.jpg")
dev.off()


ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Disease)) +
  #ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  #scale_shape_manual(values = c(4,20)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

dev.copy(jpeg,filename="PCA.jpg");
dev.off ();

# Convert to numerical expressions and save output
exprs <- exprs(eset)

ks.test(exprs, pnorm, mean(exprs), sd(exprs))

library(nortest)
ad.test(exprs)
#create Q-Q plot for both datasets
qqnorm(exprs)
qqline(exprs)
dev.copy(jpeg,filename="Normal Q-Q plot.jpg");
dev.off ();

write.csv(exprs, file = "normalized_exprs.csv")
write.table(exprs, 'normalized_exprs.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

exprs <- exprs(eset)
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
a1=gsub("\\s\\d|\\.\\s", "", c(palmieri_final@phenoData@data[["Sample number"]]))
a2=gsub("\\..*$", "", c(palmieri_final@phenoData@data[["Sample number"]]))
colnames(exprs)=as.numeric(a2)
rownames(exprs)=c(palmieri_final@featureData@data[["SYMBOL"]])
write.csv(exprs, file = "filtered_normalized_exprs.csv",row.names = TRUE, col.names = TRUE)
write.table(exprs, 'filtered_normalized_exprs.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
colnames(exprs)=palmieri_final@phenoData@data[["Time.point"]]
rownames(exprs)=c(palmieri_final@featureData@data[["SYMBOL"]])
write.csv(exprs, file = "filtered_normalized_exprs2.csv",row.names = TRUE, col.names = TRUE)
write.table(exprs, 'filtered_normalized_exprs2.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
colnames(exprs)=palmieri_final@phenoData@data[["Group.name"]]
rownames(exprs)=c(palmieri_final@featureData@data[["SYMBOL"]])
#ensemble_sym_ <- left_join(c(palmieri_final@featureData@data[["SYMBOL"]]),ensemble_sym)
rownames(expr)=c(ensemble_sym)
#exprs$SYMBOL = c(palmieri_final@featureData@data[["SYMBOL"]])
#ensemble_sym_ <- left_join(data.frame(exprs),data.frame(ensemble_sym))

write.csv(exprs, file = "filtered_normalized_exprs3.csv",row.names = TRUE, col.names = TRUE)
write.table(exprs, 'filtered_normalized_exprs3.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


# Combine two columns into one list
combined_list <- paste(palmieri_final@phenoData@data[["Group.name"]], palmieri_final@phenoData@data[["Time.point"]], sep=" ")
colnames(exprs)=combined_list
rownames(exprs)=c(palmieri_final@featureData@data[["SYMBOL"]])
write.csv(exprs, file = "filtered_normalized_exprs_Symbol_timepoint_Groupname.csv",row.names = TRUE, col.names = TRUE)
write.table(exprs, 'filtered_normalized_exprs_Symbol_timepoint_Groupname.txt', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#palmieri_final
#gene_ids <- AnnotationDbi::select(clariomdhumanhsrefseq.db, keys = rownames(palmieri_final@assayData[["exprs"]]), columns = "ENSEMBL", keytype = "PROBEID")
#aaa=merge(AA, gene_ids, by = "PROBEID", all = TRUE)
#write.csv(aaa, file = "normalized_exprs_FEMALE.csv")

# assign samples to groups and set up design matrix

#Stroke T0
#Stroke T1
#Control
#Heat-stress
palmieri_final_plot <- palmieri_final[ ,palmieri_final@phenoData@data$"Time.point"=="Control" | palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T1" | palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T0" ]
palmieri_final_young <- palmieri_final[ ,palmieri_final@phenoData@data$"Age.category"=="Adult" ]

g = palmieri_final_young@assayData[["exprs"]]
m<-mean(g)
std<-sqrt(var(g))
hist(g, density=20, breaks=20, prob=TRUE, freq = FALSE,
     xlab="Median intensitiese", 
     main="Histogram of the median intensities Young subjects")
lines(density(g), col = 4, lwd = 2)
dev.copy(jpeg,filename="Histogram_young.jpg");
dev.off ();
palmieri_final_old <- palmieri_final[ ,palmieri_final@phenoData@data$"Age.category"=="Older Adult" ]

g = palmieri_final_old@assayData[["exprs"]]
m<-mean(g)
std<-sqrt(var(g))
hist(g, density=20, breaks=20, prob=TRUE, freq = FALSE,
     xlab="Median intensitiese", 
     main="Histogram of the median intensities Old subjects")
lines(density(g), col = 4, lwd = 2)
dev.copy(jpeg,filename="Histogram_old.jpg");
dev.off ();








aggregate_strings <- function(x) {
  if (is.character(x)) {
    if (length(x) > 1) {
      return(x[1])
    } else {
      return(x)
    }
  } else {
    return(mean(x))
  }
}
palmieri_final_ <- palmieri_final[ ,palmieri_final@phenoData@data$"Time.point"=="Control" | palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T1"| palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T0"]
Subject=factor(palmieri_final_@phenoData@data$Group.number)
gs <- factor(palmieri_final_@phenoData@data$Time.point)
gender= factor(palmieri_final_@phenoData@data$Gender)
age=factor(palmieri_final_@phenoData@data$Age.category)
batch=factor(palmieri_final_@phenoData@data$Batch)
groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))
levels(gs) <- groups
palmieri_final_$group <- gs

design <- model.matrix(~gs+age+batch)
rownames(design)=palmieri_final_@phenoData@data$Time.point
colnames(design)[1:3] <- c("Stress","T0","T1")
design2<- model.matrix(~gs+age)
rownames(design2)=palmieri_final_@phenoData@data$Time.point
colnames(design2)[1:3] <- c("Stress","T0","T1")
#norm_cpn <- removeBatchEffect(palmieri_final_, batch=batch, design = design2) 
#rownames(norm_cpn)=palmieri_final_@featureData@data$SYMBOL

pheno=pData(palmieri_final_1)
mod = model.matrix(~gs+age, data=pheno)
norm_cpn <- ComBat(dat=palmieri_final_1, batch=batch, mod=mod) 
rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

exprs=exprs(palmieri_final_)
PCA_raw <- prcomp(t(exprs), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./PCA/PCA_batch_effect_PC1_PC2.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./PCA/PCA_batch_effect_PC2_PC3.jpg")
dev.off()

fit <- lmFit(norm_cpn, design2)  # fit linear model
cont.matrix <- makeContrasts(contrasts=c("T0"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)
# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_T0=top_genes_mean_table
write.csv(top_genes_mean_table, "top_genes_Stress_T0_adj_Batch.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "top_genes_Stress_T0_adj_batch.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "Shared_genes_Stress_T0_adj_bath.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B",number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Female subjets: Control- Heat Stroke T0',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="Volcano_Stress_T0_adj_bath.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=20) 

dev.copy(jpeg,filename="Barplot_Stress_T0_adj_batch.jpg");
dev.off ();
#__________________T1
cont.matrix <- makeContrasts(contrasts=c("T1"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t1=top_genes_mean_table
write.csv(top_genes_mean_table, "top_genes_Stress_T1_adj_Batch.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "top_genes_Stress_T1_adj_batch.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "Shared_genes_Stress_T1_adj_batch.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B", number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Female subjets: Control- Heat Stroke T1',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="Volcano_Stress_T1_adj_batch.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="Barplot_Stress_T1_adj_batch.jpg");
dev.off ();

a <- list('DE T0' = top_genes_mean_table_T0$SYMBOL,
          'DE T1' = top_genes_mean_table_t1$SYMBOL)
venn <- ggvenn(a)
# Add a title to the Venn diagram
venn <- venn + ggtitle("DE Genes")
# Print the Venn diagram
print(venn)
dev.copy(jpeg,filename="DE_genes_T0_T1.jpg");
dev.off ()




# Create the input vectors.

colors = c("green", "red" )

stnames <-  c("HS-T0", "HS-T1") 

subjects <-  c("Upregulated", "Downregulated")       

# Create the matrix of the values.

Values <- matrix( c(2529,2054,   3003 ,2422 ) , nrow=2, ncol=2, byrow=TRUE)  

# Create the bar chart

barplot(Values, names.arg = stnames, ylab = "Number of DE Genes", col = colors, beside=TRUE ) 

#Ad an optional legend

legend("bottomright", subjects, cex = 0.7, fill = colors) 

dev.copy(jpeg,filename="DE_genes_T0_T1_barplot.jpg");
dev.off ()


aggregate_strings <- function(x) {
  if (is.character(x)) {
    if (length(x) > 1) {
      return(x[1])
    } else {
      return(x)
    }
  } else {
    return(mean(x))
  }
}
palmieri_final_ <- palmieri_final[ ,palmieri_final@phenoData@data$"Time.point"=="Control" | palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T1"| palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T0"]
palmieri_final_1 <- palmieri_final_[ ,palmieri_final_@phenoData@data$"Gender"=="Female"]

Subject=factor(palmieri_final_1@phenoData@data$Group.number)
gs <- factor(palmieri_final_1@phenoData@data$Time.point)
gender= factor(palmieri_final_1@phenoData@data$Gender)
age=factor(palmieri_final_1@phenoData@data$Age.category)
batch=factor(palmieri_final_1@phenoData@data$Batch)
label=factor(palmieri_final_1@phenoData@data[["Group.name"]])
label=paste(palmieri_final@phenoData@data[["Group.name"]],palmieri_final@phenoData@data[["Time.point"]],palmieri_final@phenoData@data[["Age.category"]])

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))
levels(gs) <- groups
palmieri_final_1$group <- gs

design <- model.matrix(~gs+age+batch)
colnames(design)[1:3] <- c("Stress","T0","T1")
rownames(design)=palmieri_final_1@phenoData@data$Time.point
#covariate <- model.matrix(~age)
#batch = model.matrix(~age+batch)
design2<- model.matrix(~gs+age)
colnames(design2)[1:3] <- c("Stress","T0","T1")
rownames(design2)=palmieri_final_1@phenoData@data$Time.point

pheno=pData(palmieri_final_1)
mod = model.matrix(~gs+age, data=pheno)
norm_cpn <- ComBat(dat=palmieri_final_1, batch=batch, mod=mod) 
rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL
#norm_cpn <- removeBatchEffect(palmieri_final_1, batch=batch, design = design2) 
#rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

exprs=exprs(palmieri_final_1)
PCA_raw <- prcomp(t(exprs), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./FEMALE/PCA_batch_PC1_PC2.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./FEMALE/PCA_batch_PC2_PC3.jpg")
dev.off()


PCA_raw <- prcomp(t(norm_cpn), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./FEMALE/PCA_batch_PC1_PC2_remove_batch_eff_Combat.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./FEMALE/PCA_PC2_PC3_remove_batch_eff_Combat.jpg")
dev.off()

fit <- lmFit(norm_cpn, design2)  # fit linear model
cont.matrix <- makeContrasts(contrasts=c("T0"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_T0_female=top_genes_mean_table
write.csv(top_genes_mean_table, "./FEMALE/top_genes_Stress_T0_FEMALE_adj_Batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "./FEMALE/top_genes_Stress_T0_FEMALE_adj_batch_Combat.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./FEMALE/Shared_genes_Stress_T0_FEMALE_adj_bath_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B",number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Female subjets: Control- Heat Stroke T0',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./FEMALE/Volcano_Stress_T0_FEMALE_adj_bath_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./FEMALE/Barplot_Stress_T0_FEMALE_adj_batch_Combat.jpg");
dev.off ();


#___________________T1

cont.matrix <- makeContrasts(contrasts=c("T1"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t1_female=top_genes_mean_table
write.csv(top_genes_mean_table, "./FEMALE/top_genes_Stress_T1_FEMALE_adj_Batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "./FEMALE/top_genes_Stress_T1_FEMALE_adj_batch_Combat.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./FEMALE/Shared_genes_Stress_T1_FEMALE_adj_batch_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B", number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Female subjets: Control- Heat Stroke T1',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./FEMALE/Volcano_Stress_T1_FEMALE_adj_batch_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./FEMALE/Barplot_Stress_T1_FEMALE_adj_batch_Combat.jpg");
dev.off ();

#______________MALE

palmieri_final_ <- palmieri_final[ ,palmieri_final@phenoData@data$"Time.point"=="Control" | palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T1"| palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T0"]
palmieri_final_1 <- palmieri_final_[ ,palmieri_final_@phenoData@data$"Gender"=="Male"]

Subject=factor(palmieri_final_1@phenoData@data$Group.number)
gs <- factor(palmieri_final_1@phenoData@data$Time.point)
gender= factor(palmieri_final_1@phenoData@data$Gender)
age=factor(palmieri_final_1@phenoData@data$Age.category)
batch=factor(palmieri_final_1@phenoData@data$Batch)

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))
levels(gs) <- groups
palmieri_final_1$group <- gs

design <- model.matrix(~gs+age+batch)
colnames(design)[1:3] <- c("Stress","T0","T1")
rownames(design)=palmieri_final_1@phenoData@data$Time.point

design2<- model.matrix(~gs+age)
colnames(design2)[1:3] <- c("Stress","T0","T1")
rownames(design2)=palmieri_final_1@phenoData@data$Time.point
#norm_cpn <- removeBatchEffect(palmieri_final_1, batch=batch, design = design2) 
#rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

pheno=pData(palmieri_final_1)
mod = model.matrix(~gs+age, data=pheno)
norm_cpn <- ComBat(dat=palmieri_final_1, batch=batch, mod=mod) 
rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

exprs=exprs(palmieri_final_1)
PCA_raw <- prcomp(t(exprs), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./MALE/PCA_batch_PC1_PC2.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./MALE/PCA_batch_PC2_PC3.jpg")
dev.off()


PCA_raw <- prcomp(t(norm_cpn), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./MALE/PCA_batch_PC1_PC2_remove_batch_eff_Combat.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./MALE/PCA_PC2_PC3_remove_batch_eff_Combat.jpg")
dev.off()

fit <- lmFit(norm_cpn, design2)  # fit linear model

cont.matrix <- makeContrasts(contrasts="T0", levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", number=Inf)
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05, number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t0_male = top_genes_mean_table
write.csv(top_genes_mean_table, "./MALE/top_genes_stress_T0_MALE_adj_batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "./MALE/top_genes_Stress_T0_MALE_adj_batch_Combat.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./MALE/Shared_genes_stress_T0_MALE_adj_batch_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B", number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Male subjets: Control- Heat Stroke T0',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./MALE/Volcano_stress_T0_MALE_adj_batch_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./MALE/Barplot_stress_T0_MALE_adj_batch_Combat.jpg");
dev.off ();

#_______________________T1

cont.matrix <- makeContrasts(contrasts="T1", levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t1_male=top_genes_mean_table
write.csv(top_genes_mean_table, "./MALE/top_genes_Stress_T1_MALE_adj_batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes, gene_ids, by = "SYMBOL", all = TRUE)
Ref <- read_excel("Hest_stroke_genes.xlsx")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./MALE/Shared_genes_Stress_T1_MALE_adj_batch_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B",number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Male subjets: Control- Heat Stroke T1',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./MALE/Volcano_Stress_T1_MALE_adj_batch_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./MALE/Barplot_Stress_T1_MALE_adj_batch_Combat.jpg");
dev.off ();
a=1

#________________________YOUNG

palmieri_final_ <- palmieri_final[ ,palmieri_final@phenoData@data$"Time.point"=="Control" | palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T1"| palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T0"]
palmieri_final_1 <- palmieri_final_[ ,palmieri_final_@phenoData@data$Age.category=="Young"]

Subject=factor(palmieri_final_1@phenoData@data$Group.number)
gs <- factor(palmieri_final_1@phenoData@data$Time.point)
gender= factor(palmieri_final_1@phenoData@data$Gender)
age=factor(palmieri_final_1@phenoData@data$Age.category)
batch=factor(palmieri_final_1@phenoData@data$Batch)
label=factor(palmieri_final_1@phenoData@data[["Group.name"]])
label=paste(palmieri_final@phenoData@data[["Group.name"]],palmieri_final@phenoData@data[["Time.point"]],palmieri_final@phenoData@data[["Age.category"]])

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))
levels(gs) <- groups
palmieri_final_1$group <- gs

design <- model.matrix(~gs+gender+batch)
colnames(design)[1:3] <- c("Stress","T0","T1")
rownames(design)=palmieri_final_1@phenoData@data$Time.point
#covariate <- model.matrix(~age)
#batch = model.matrix(~age+batch)
design2<- model.matrix(~gs+gender)
rownames(design2)=palmieri_final_1@phenoData@data$Time.point
colnames(design2)[1:3] <- c("Stress","T0","T1")
#norm_cpn <- removeBatchEffect(palmieri_final_1, batch=batch, design = design2) 
#rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

pheno=pData(palmieri_final_1)
mod = model.matrix(~gs+gender, data=pheno)
norm_cpn <- ComBat(dat=palmieri_final_1, batch=batch, mod=mod) 
rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

exprs=exprs(palmieri_final_1)
PCA_raw <- prcomp(t(exprs), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./YOUNG/PCA_batch_PC1_PC2.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./YOUNG/PCA_batch_PC2_PC3.jpg")
dev.off()


PCA_raw <- prcomp(t(norm_cpn), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./YOUNG/PCA_batch_PC1_PC2_remove_batch_eff_Combat.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./YOUNG/PCA_PC2_PC3_remove_batch_eff_Combat.jpg")
dev.off()

fit <- lmFit(norm_cpn, design2)  # fit linear model
cont.matrix <- makeContrasts(contrasts=c("T0"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t0_young=top_genes_mean_table
write.csv(top_genes_mean_table, "./YOUNG/top_genes_Stress_T0_YOUNG_adj_Batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "./YOUNG/top_genes_Stress_T0_YOUNG_adj_batch_Combat.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./YOUNG/Shared_genes_Stress_T0_YOUNG_adj_bath_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B",number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Young subjets: Control- Heat Stroke T0',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./YOUNG/Volcano_Stress_T0_YOUNG_adj_bath_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./YOUNG/Barplot_Stress_T0_YOUNG_adj_batch_Combat.jpg");
dev.off ();


#___________________T1

cont.matrix <- makeContrasts(contrasts=c("T1"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t1_young=top_genes_mean_table
write.csv(top_genes_mean_table, "./YOUNG/top_genes_Stress_T1_YOUNG_adj_Batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "./YOUNG/top_genes_Stress_T1_YOUNG_adj_batch_Combat.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./YOUNG/Shared_genes_Stress_T1_YOUNG_adj_batch_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B", number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Young subjets: Control- Heat Stroke T1',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./YOUNG/Volcano_Stress_T1_YOUNG_adj_batch_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./YOUNG/Barplot_Stress_T1_YOUNG_adj_batch_Combat.jpg");
dev.off ();

#______________________OLD

palmieri_final_ <- palmieri_final[ ,palmieri_final@phenoData@data$"Time.point"=="Control" | palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T1"| palmieri_final@phenoData@data$"Time.point"=="Heat Stroke T0"]
palmieri_final_1 <- palmieri_final_[ ,palmieri_final_@phenoData@data$Age.category=="Old"]

Subject=factor(palmieri_final_1@phenoData@data$Group.number)
gs <- factor(palmieri_final_1@phenoData@data$Time.point)
gender= factor(palmieri_final_1@phenoData@data$Gender)
age=factor(palmieri_final_1@phenoData@data$Age.category)
batch=factor(palmieri_final_1@phenoData@data$Batch)
label=factor(palmieri_final_1@phenoData@data[["Group.name"]])
label=paste(palmieri_final@phenoData@data[["Group.name"]],palmieri_final@phenoData@data[["Time.point"]],palmieri_final@phenoData@data[["Age.category"]])

groups <-make.names(c("Stress","T1","T0"))# make.names(c("control","control","T0","T0"))
levels(gs) <- groups
palmieri_final_1$group <- gs

design <- model.matrix(~gs+gender+batch)
colnames(design)[1:3] <- c("Stress","T0","T1")
#covariate <- model.matrix(~age)
#batch = model.matrix(~age+batch)
design2<- model.matrix(~gs+gender)
colnames(design2)[1:3] <- c("Stress","T0","T1")
#norm_cpn <- removeBatchEffect(palmieri_final_1, batch=batch, design = design2) 
#rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

pheno=pData(palmieri_final_1)
mod = model.matrix(~gs+gender, data=pheno)
norm_cpn <- ComBat(dat=palmieri_final_1, batch=batch, mod=mod) 
rownames(norm_cpn)=palmieri_final_1@featureData@data$SYMBOL

exprs=exprs(palmieri_final_1)
PCA_raw <- prcomp(t(exprs), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./OLD/PCA_batch_PC1_PC2.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./OLD/PCA_batch_PC2_PC3.jpg")
dev.off()


PCA_raw <- prcomp(t(norm_cpn), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],
                     Batch = batch,
                     Disease = gs)

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC1, PC2))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC1, y = PC2, color = Batch) +  # Color by Disease
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +  # Assuming percentVar is defined
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +  # Assuming percentVar is defined
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)# Assuming sd_ratio is defined

dev.copy(jpeg, filename = "./OLD/PCA_batch_PC1_PC2_remove_batch_eff_Combat.jpg")
dev.off()

dataGG_clean <- dataGG %>%
  drop_na() %>%
  group_by(Batch) %>%
  dplyr::slice(chull(PC2, PC3))  # Calculate the convex hull for each group

ggplot(dataGG) +
  aes(x = PC2, y = PC3, color = Batch) +
  geom_point(aes(shape = Disease)) +
  geom_polygon(data = dataGG_clean,
               aes(fill = Disease, color = NULL),  # Remove color mapping from geom_polygon
               alpha = 0.3,
               show.legend = FALSE) +
  xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ylab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio)

dev.copy(jpeg, filename = "./OLD/PCA_PC2_PC3_remove_batch_eff_Combat.jpg")
dev.off()

fit <- lmFit(norm_cpn, design2)  # fit linear model
cont.matrix <- makeContrasts(contrasts=c("T0"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t0_old=top_genes_mean_table
write.csv(top_genes_mean_table, "./OLD/top_genes_Stress_T0_OLD_adj_Batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "./OLD/top_genes_Stress_T0_OLD_adj_batch_Combat.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./OLD/Shared_genes_Stress_T0_OLD_adj_bath_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B",number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Old subjets: Control- Heat Stroke T0',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./OLD/Volcano_Stress_T0_OLD_adj_bath_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./OLD/Barplot_Stress_T0_OLD_adj_batch_Combat.jpg");
dev.off ();


#___________________T1

cont.matrix <- makeContrasts(contrasts=c("T1"), levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
DE_genes2 <- decideTests(fit2)
summary(DE_genes2)

# Get the top 10 deferentially expressed genes
top_genes <- topTable(fit2, adjust="BH", sort.by="B", p.value=0.05,number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)
top_genes_mean_table_t1_old=top_genes_mean_table
write.csv(top_genes_mean_table, "./OLD/top_genes_Stress_T1_OLD_adj_Batch_Combat.csv")

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_mean_table$SYMBOL, columns = "ENTREZID", keytype = "SYMBOL")
gene_list <- merge(top_genes_mean_table, gene_ids, by = "SYMBOL", all = TRUE)

Ref <- read_excel("Hest_stroke_genes.xlsx")
write.csv(gene_list, "./OLD/top_genes_Stress_T1_OLD_adj_batch_Combat.csv")
int=list(intersect(Ref$SYMBOL, top_genes_mean_table$SYMBOL))
merged_table <- data.frame(merge(Ref, top_genes_mean_table, by = "SYMBOL"))
write.csv(merged_table, "./OLD/Shared_genes_Stress_T1_OLD_adj_batch_Combat.csv")

top_genes <- topTable(fit2, adjust="BH", sort.by="B", number=Inf)
top_genes=na.omit(top_genes)
colnames(top_genes)[1]="SYMBOL"
top_genes_mean_table <- aggregate(. ~ SYMBOL, data = top_genes, FUN = aggregate_strings)

merged_table$logFC=as.numeric(unlist(merged_table$logFC))
merged_table$adj.P.Val=as.numeric(unlist(merged_table$adj.P.Val))
top_genes_mean_table$logFC=as.numeric(unlist(top_genes_mean_table$logFC))
top_genes_mean_table$adj.P.Val=as.numeric(unlist(top_genes_mean_table$adj.P.Val))

EnhancedVolcano(top_genes_mean_table,
                lab = top_genes_mean_table$SYMBOL,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Old subjets: Control- Heat Stroke T1',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0)
dev.copy(jpeg,filename="./OLD/Volcano_Stress_T1_OLD_adj_batch_Combat.jpg");
dev.off ();

OrgDb='org.Hs.eg.db'
ego <- enrichGO(gene_ids$ENTREZID, OrgDb, ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
#Bar plot of enriched terms.s
barplot(ego, showCategory=15) 

dev.copy(jpeg,filename="./OLD/Barplot_Stress_T1_OLD_adj_batch_Combat.jpg");
dev.off ();

top_genes_mean_table_t0_female_up <- top_genes_mean_table_T0_female[top_genes_mean_table_T0_female$logFC > 0, ]
top_genes_mean_table_t0_female_up = top_genes_mean_table_t0_female_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t0_female_down <- top_genes_mean_table_T0_female[top_genes_mean_table_T0_female$logFC < 0, ]
top_genes_mean_table_t0_female_down = top_genes_mean_table_t0_female_down[, c('SYMBOL',"logFC","adj.P.Val")]

top_genes_mean_table_t0_male_up <- top_genes_mean_table_t0_male[top_genes_mean_table_t0_male$logFC > 0, ]
top_genes_mean_table_t0_male_up = top_genes_mean_table_t0_male_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t0_male_down <- top_genes_mean_table_t0_male[top_genes_mean_table_t0_male$logFC < 0, ]
top_genes_mean_table_t0_male_down = top_genes_mean_table_t0_male_down[, c('SYMBOL',"logFC","adj.P.Val")]

top_genes_mean_table_t0_young_up <- top_genes_mean_table_t0_young[top_genes_mean_table_t0_young$logFC > 0, ]
top_genes_mean_table_t0_young_up = top_genes_mean_table_t0_young_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t0_young_down <- top_genes_mean_table_t0_young[top_genes_mean_table_t0_young$logFC < 0, ]
top_genes_mean_table_t0_young_down = top_genes_mean_table_t0_young_down[, c('SYMBOL',"logFC","adj.P.Val")]

top_genes_mean_table_t0_old_up <- top_genes_mean_table_t0_old[top_genes_mean_table_t0_old$logFC > 0, ]
top_genes_mean_table_t0_old_up = top_genes_mean_table_t0_old_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t0_old_down <- top_genes_mean_table_t0_old[top_genes_mean_table_t0_old$logFC < 0, ]
top_genes_mean_table_t0_old_down = top_genes_mean_table_t0_old_down[, c('SYMBOL',"logFC","adj.P.Val")]


list1 <- top_genes_mean_table_t0_female_up$SYMBOL
list2 <- top_genes_mean_table_t0_female_down$SYMBOL
list3 <- top_genes_mean_table_t0_male_up$SYMBOL
list4 <- top_genes_mean_table_t0_male_down$SYMBOL
list5 <- top_genes_mean_table_t0_young_up$SYMBOL
list6 <- top_genes_mean_table_t0_young_down$SYMBOL
list7 <- top_genes_mean_table_t0_old_up$SYMBOL
list8 <- top_genes_mean_table_t0_old_down$SYMBOL

# Determine the maximum length among the lists
max_length <- max(length(list1), length(list2), length(list3), length(list4), length(list5), length(list6), length(list7), length(list8))

# Pad the shorter lists with NA values to match the maximum length
list1 <- c(list1, rep(NA, max_length - length(list1)))
list2 <- c(list2, rep(NA, max_length - length(list2)))
list3 <- c(list3, rep(NA, max_length - length(list3)))
list4 <- c(list4, rep(NA, max_length - length(list4)))
list5 <- c(list5, rep(NA, max_length - length(list5)))
list6 <- c(list6, rep(NA, max_length - length(list6)))
list7 <- c(list7, rep(NA, max_length - length(list7)))
list8 <- c(list8, rep(NA, max_length - length(list8)))

# Combine the lists into a data frame
result_table_sym <- data.frame(list1, list2, list3, list4, list5, list6, list7, list8)
colnames(result_table_sym) <- c("Female Up", "Female Down", "Male Up", "Male Down", "Young Up", "Young Down", "Old Up", "Old Down")

list1 <- top_genes_mean_table_t0_female_up$logFC
list2 <- top_genes_mean_table_t0_female_down$logFC
list3 <- top_genes_mean_table_t0_male_up$logFC
list4 <- top_genes_mean_table_t0_male_down$logFC
list5 <- top_genes_mean_table_t0_young_up$logFC
list6 <- top_genes_mean_table_t0_young_down$logFC
list7 <- top_genes_mean_table_t0_old_up$logFC
list8 <- top_genes_mean_table_t0_old_down$logFC

# Pad the shorter lists with NA values to match the maximum length
list1 <- c(list1, rep(NA, max_length - length(list1)))
list2 <- c(list2, rep(NA, max_length - length(list2)))
list3 <- c(list3, rep(NA, max_length - length(list3)))
list4 <- c(list4, rep(NA, max_length - length(list4)))
list5 <- c(list5, rep(NA, max_length - length(list5)))
list6 <- c(list6, rep(NA, max_length - length(list6)))
list7 <- c(list7, rep(NA, max_length - length(list7)))
list8 <- c(list8, rep(NA, max_length - length(list8)))

# Combine the lists into a data frame
result_table_logFC <- data.frame(list1, list2, list3, list4, list5, list6, list7, list8)
colnames(result_table_logFC) <- c("Female Up", "Female Down", "Male Up", "Male Down", "Young Up", "Young Down", "Old Up", "Old Down")

list1 <- top_genes_mean_table_t0_female_up$adj.P.Val
list2 <- top_genes_mean_table_t0_female_down$adj.P.Val
list3 <- top_genes_mean_table_t0_male_up$adj.P.Val
list4 <- top_genes_mean_table_t0_male_down$adj.P.Val
list5 <- top_genes_mean_table_t0_young_up$adj.P.Val
list6 <- top_genes_mean_table_t0_young_down$adj.P.Val
list7 <- top_genes_mean_table_t0_old_up$adj.P.Val
list8 <- top_genes_mean_table_t0_old_down$adj.P.Val

# Pad the shorter lists with NA values to match the maximum length
list1 <- c(list1, rep(NA, max_length - length(list1)))
list2 <- c(list2, rep(NA, max_length - length(list2)))
list3 <- c(list3, rep(NA, max_length - length(list3)))
list4 <- c(list4, rep(NA, max_length - length(list4)))
list5 <- c(list5, rep(NA, max_length - length(list5)))
list6 <- c(list6, rep(NA, max_length - length(list6)))
list7 <- c(list7, rep(NA, max_length - length(list7)))
list8 <- c(list8, rep(NA, max_length - length(list8)))

# Combine the lists into a data frame
result_table_padj <- data.frame(list1, list2, list3, list4, list5, list6, list7, list8)
colnames(result_table_padj) <- c("Female Up", "Female Down", "Male Up", "Male Down", "Young Up", "Young Down", "Old Up", "Old Down")
#write.csv(result_table_padj, "result_table_T0.csv")
output_file <- 'Up_down_categories_Combat.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Genes')
writeData(wb, sheet = 'Genes', x = result_table_sym)
addWorksheet(wb, sheetName = 'LogFC')
writeData(wb, sheet = 'LogFC', x = result_table_logFC)
addWorksheet(wb, sheetName = 'p_adj')
writeData(wb, sheet = 'p_adj', x = result_table_padj)
saveWorkbook(wb, output_file)


top_genes_mean_table_t1_female_up <- top_genes_mean_table_t1_female[top_genes_mean_table_t1_female$logFC > 0, ]
top_genes_mean_table_t1_female_up = top_genes_mean_table_t1_female_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t1_female_down <- top_genes_mean_table_t1_female[top_genes_mean_table_t1_female$logFC < 0, ]
top_genes_mean_table_t1_female_down = top_genes_mean_table_t1_female_down[, c('SYMBOL',"logFC","adj.P.Val")]

top_genes_mean_table_t1_male_up <- top_genes_mean_table_t1_male[top_genes_mean_table_t1_male$logFC > 0, ]
top_genes_mean_table_t1_male_up = top_genes_mean_table_t1_male_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t1_male_down <- top_genes_mean_table_t1_male[top_genes_mean_table_t1_male$logFC < 0, ]
top_genes_mean_table_t1_male_down = top_genes_mean_table_t1_male_down[, c('SYMBOL',"logFC","adj.P.Val")]

top_genes_mean_table_t1_young_up <- top_genes_mean_table_t1_young[top_genes_mean_table_t1_young$logFC > 0, ]
top_genes_mean_table_t1_young_up = top_genes_mean_table_t1_young_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t1_young_down <- top_genes_mean_table_t1_young[top_genes_mean_table_t1_young$logFC < 0, ]
top_genes_mean_table_t1_young_down = top_genes_mean_table_t1_young_down[, c('SYMBOL',"logFC","adj.P.Val")]

top_genes_mean_table_t1_old_up <- top_genes_mean_table_t1_old[top_genes_mean_table_t1_old$logFC > 0, ]
top_genes_mean_table_t1_old_up = top_genes_mean_table_t1_old_up[, c('SYMBOL',"logFC","adj.P.Val")]
top_genes_mean_table_t1_old_down <- top_genes_mean_table_t1_old[top_genes_mean_table_t1_old$logFC < 0, ]
top_genes_mean_table_t1_old_down = top_genes_mean_table_t1_old_down[, c('SYMBOL',"logFC","adj.P.Val")]

list1 <- top_genes_mean_table_t1_female_up$SYMBOL
list2 <- top_genes_mean_table_t1_female_down$SYMBOL
list3 <- top_genes_mean_table_t1_male_up$SYMBOL
list4 <- top_genes_mean_table_t1_male_down$SYMBOL
list5 <- top_genes_mean_table_t1_young_up$SYMBOL
list6 <- top_genes_mean_table_t1_young_down$SYMBOL
list7 <- top_genes_mean_table_t1_old_up$SYMBOL
list8 <- top_genes_mean_table_t1_old_down$SYMBOL

# Determine the maximum length among the lists
max_length <- max(length(list1), length(list2), length(list3), length(list4), length(list5), length(list6), length(list7), length(list8))

# Pad the shorter lists with NA values to match the maximum length
list1 <- c(list1, rep(NA, max_length - length(list1)))
list2 <- c(list2, rep(NA, max_length - length(list2)))
list3 <- c(list3, rep(NA, max_length - length(list3)))
list4 <- c(list4, rep(NA, max_length - length(list4)))
list5 <- c(list5, rep(NA, max_length - length(list5)))
list6 <- c(list6, rep(NA, max_length - length(list6)))
list7 <- c(list7, rep(NA, max_length - length(list7)))
list8 <- c(list8, rep(NA, max_length - length(list8)))

# Combine the lists into a data frame
result_table_sym <- data.frame(list1, list2, list3, list4, list5, list6, list7, list8)
colnames(result_table_sym) <- c("Female Up", "Female Down", "Male Up", "Male Down", "Young Up", "Young Down", "Old Up", "Old Down")

list1 <- top_genes_mean_table_t1_female_up$logFC
list2 <- top_genes_mean_table_t1_female_down$logFC
list3 <- top_genes_mean_table_t1_male_up$logFC
list4 <- top_genes_mean_table_t1_male_down$logFC
list5 <- top_genes_mean_table_t1_young_up$logFC
list6 <- top_genes_mean_table_t1_young_down$logFC
list7 <- top_genes_mean_table_t1_old_up$logFC
list8 <- top_genes_mean_table_t1_old_down$logFC

# Pad the shorter lists with NA values to match the maximum length
list1 <- c(list1, rep(NA, max_length - length(list1)))
list2 <- c(list2, rep(NA, max_length - length(list2)))
list3 <- c(list3, rep(NA, max_length - length(list3)))
list4 <- c(list4, rep(NA, max_length - length(list4)))

list5 <- c(list5, rep(NA, max_length - length(list5)))
list6 <- c(list6, rep(NA, max_length - length(list6)))
list7 <- c(list7, rep(NA, max_length - length(list7)))
list8 <- c(list8, rep(NA, max_length - length(list8)))

# Combine the lists into a data frame
result_table_logFC <- data.frame(list1, list2, list3, list4, list5, list6, list7, list8)
colnames(result_table_logFC) <- c("Female Up", "Female Down", "Male Up", "Male Down", "Young Up", "Young Down", "Old Up", "Old Down")

list1 <- top_genes_mean_table_t1_female_up$adj.P.Val
list2 <- top_genes_mean_table_t1_female_down$adj.P.Val
list3 <- top_genes_mean_table_t1_male_up$adj.P.Val
list4 <- top_genes_mean_table_t1_male_down$adj.P.Val
list5 <- top_genes_mean_table_t1_young_up$adj.P.Val
list6 <- top_genes_mean_table_t1_young_down$adj.P.Val
list7 <- top_genes_mean_table_t1_old_up$adj.P.Val
list8 <- top_genes_mean_table_t1_old_down$adj.P.Val

# Pad the shorter lists with NA values to match the maximum length
list1 <- c(list1, rep(NA, max_length - length(list1)))
list2 <- c(list2, rep(NA, max_length - length(list2)))
list3 <- c(list3, rep(NA, max_length - length(list3)))
list4 <- c(list4, rep(NA, max_length - length(list4)))
list5 <- c(list5, rep(NA, max_length - length(list5)))
list6 <- c(list6, rep(NA, max_length - length(list6)))
list7 <- c(list7, rep(NA, max_length - length(list7)))
list8 <- c(list8, rep(NA, max_length - length(list8)))

# Combine the lists into a data frame
result_table_padj <- data.frame(list1, list2, list3, list4, list5, list6, list7, list8)
colnames(result_table_padj) <- c("Female Up", "Female Down", "Male Up", "Male Down", "Young Up", "Young Down", "Old Up", "Old Down")

output_file <- 'Up_down_categories_T1_Combat.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Genes')
writeData(wb, sheet = 'Genes', x = result_table_sym)
addWorksheet(wb, sheetName = 'LogFC')
writeData(wb, sheet = 'LogFC', x = result_table_logFC)
addWorksheet(wb, sheetName = 'p_adj')
writeData(wb, sheet = 'p_adj', x = result_table_padj)
saveWorkbook(wb, output_file)
















# Load necessary libraries
library(readxl)
library(dplyr)
library(gplots)
library(ggplot2)
library(GEOquery)
library(clusterProfiler)
library(clariomdhumanhsrefseq.db)
library(clariomdhumanhsrefseqcdf)
library(pd.clariomdhuman.hs.refseq)
library(org.Hs.eg.db) # Added library for AnnotationDbi::select

# Function to read and process data
read_proteomics_data <- function(file_path) {
  proteomics_data <- read.csv(file_path)
  data_ <- AnnotationDbi::select(org.Hs.eg.db, keys = proteomics_data$X, columns = "SYMBOL", keytype = "UNIPROT")
  colnames(data_)[1] <- "X"
  proteomics_data <- left_join(data_, proteomics_data, by = "X")
  proteomics_data <- proteomics_data[!duplicated(proteomics_data$SYMBOL), ]
  rownames(proteomics_data) <- proteomics_data$SYMBOL
  proteomics_data <- proteomics_data %>% select(-X, -SYMBOL)
  return(proteomics_data)
}

# Example file paths
proteomics_file <- "data/normalized_data_Standard_Deviation.csv"

# Read data
transcriptomic_data <- data.frame(exprs)
transcriptomic_data <- transcriptomic_data[, !grepl("CD", names(transcriptomic_data))]
proteomics_data <- read_proteomics_data(proteomics_file)

# Find common genes and align data
common_genes <- intersect(proteomics_data$SYMBOL, rownames(transcriptomic_data))
trans <- transcriptomic_data[common_genes, ]
prot <- proteomics_data[common_genes, ]

# Remove duplicates
prot <- prot[!duplicated(prot$SYMBOL), ]
rownames(prot) <- prot$SYMBOL

# Reorder columns
trans <- trans[colnames(prot)]

# Perform correlation tests
correlation_results <- lapply(1:ncol(prot), function(i) {
  cor.test(prot[, i], trans[, i])
})

# Write correlation results to file
write.csv(correlation_results, "results/correlation_all_genes.csv")

# Subset and correlate specific time points
subset_timepoints <- function(trans, prot, timepoint) {
  trans_subset <- trans[, grepl(timepoint, names(trans))]
  prot_subset <- prot[, grepl(timepoint, names(prot))]
  return(cor(t(prot_subset), t(trans_subset)))
}

correlation_t0_t0 <- subset_timepoints(trans, prot, "T0")
correlation_t1_t1 <- subset_timepoints(trans, prot, "T1")

# Write specific correlation results
write.csv(correlation_t0_t0, "results/correlation_t0_t0.csv")
write.csv(correlation_t1_t1, "results/correlation_t1_t1.csv")

# Differential expression analysis
process_DE_data <- function(file_path) {
  de_data <- read.csv(file_path)
  de_data <- de_data[!duplicated(de_data$SYMBOL), ]
  rownames(de_data) <- de_data$SYMBOL
  return(de_data %>% select(-SYMBOL))
}

proteomic_DE_file <- "data/proteomic_DE.csv"
transcriptomic_DE_file <- "data/transcriptomic_DE.csv"

proteomic_DE <- process_DE_data(proteomic_DE_file)
transcriptomic_DE <- process_DE_data(transcriptomic_DE_file)

comp <- intersect(rownames(proteomic_DE), rownames(transcriptomic_DE))

# Correlation of differential expression data
correlation_DE <- cor(transcriptomic_DE[comp, ], proteomic_DE[comp, ])

# Save correlation results
write.csv(correlation_DE, "results/correlation_DE.csv")

library(RnBeads)
library("RnBeads.hg38")
library(RnBeads)
library(doParallel)
library("RnBeads.hg38")
library(shiny)
library(shinyjs)
library(plotrix)
library("ggbio")
library("missMethyl")
library("GOstats")

setwd("/Users/mariacastillo/Desktop")
data.dir <- "/Users/mariacastillo/Desktop/RnBeads"
sample.annotation <- file.path(data.dir, "file_list.csv")
analysis.dir <- "/Users/mariacastillo/Desktop/Report_13_18_23"
report.dir <- file.path(analysis.dir, "reports")
rnb.options("assembly" = "hg38")
rnb.options("filtering.greedycut"=FALSE)
rnb.options("exploratory.intersample"=FALSE)
rnb.options("distribution.subsample"=100000)
rnb.options("exploratory.region.profiles"=NULL)
rnb.options("qc.coverage.plots"=FALSE)
rnb.options("inference.age.prediction"=FALSE)
rnb.options("identifiers.column"="sampleID")
rnb.options("import.bed.style" = "bismarkCov")
rnb.options("normalization" = TRUE)
rnb.options("filtering.low.coverage.masking" = FALSE)
rnb.options("filtering.high.coverage.outliers" = FALSE)
rnb.options("filtering.sex.chromosomes.removal"=TRUE)
rnb.options("inference" = FALSE)
rnb.options("inference.age.column" = "Age")
rnb.options("exploratory.columns" = c("differentiation_level_1","Predicted Sex","Age"))
rnb.options("differential.variability" = FALSE)
rnb.options("differential.variability.method" = "diffVar")
rnb.options("differential.comparison.columns"=c("differentiation_level_1","Predicted Sex","Age") )
rnb.options("differential.enrichment.go"=TRUE)
rnb.options("differential.report.sites" = TRUE)
rnb.options("export.to.bed" = TRUE)
rnb.options("export.to.csv"=TRUE)
export.types = c("sites","tiling","genes","promoters","cpgislands")
 
#result <- rnb.run.import(data.source=files, data.type="bs.bed.dir", dir.reports=report.dir)

rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation, data.dir=data.dir, data.type="bs.bed.dir")
#?rnb.run.analysis
load.rnb.set("/Users/mariacastillo/Desktop/Methylation_Results/RnBeads_reports/rnbSet_preprocessed", temp.dir = tempdir())


library(RnBeads)
library(wordcloud)
rnb.run.dj()
rnb.set=load.rnb.set("/Users/mariacastillo/Desktop/Methylation_Results/RnBeads_reports/rnbSet_preprocessed", temp.dir = tempdir())
diffmeth=load.rnb.diffmeth("/Users/mariacastillo/Desktop/HEATSTROKE/Methylation/results_final/reports_ALL/differential_methylation_data/differential_rnbDiffMeth")
rnb.options("differential.variability"=TRUE)
rnb.options("differential.comparison.columns"=c("differentiation_level_1"), "columns.pairing"=c("differentiation_level_1"="Predicted Sex"))
rnb.options("differential.enrichment.go"=TRUE)
rnb.options("differential.enrichment.lola"=FALSE)
rnb.options(export.to.csv=TRUE)
rnb.run.differential(rnb.set, report.dir)
enrich.go <- performGoEnrichment.diffMeth(rnb.set, diffmeth, verbose=TRUE)
enrich.table.go <- enrich.go[["region"]][["Heatstress_control vs. Heatstroke_patient (based on differentiation_level_1)"]][["BP"]][["promoters"]][["rankCut_500"]][["hypo"]]
class(enrich.table.go)
summary(enrich.table.go)
genes_enrich_go=data.frame(enrich.go[["region"]][["Heatstress_control vs. Heatstroke_patient (based on differentiation_level_1)"]][["BP"]][["genes"]][["rankCut_autoSelect"]][["hyper"]]@geneIds)
library("writexl")
write_xlsx(genes_enrich_go,"genes_enrich_go.xlsx")


cmp.cols <- "Sample_Group"
reg.types <- c("genes", "promoters")
diffmeth <- rnb.execute.computeDiffMeth(rnb.set, cmp.cols, region.types=reg.types)

enrich.go=load.rnb.set("/Users/mariacastillo/Desktop/Report_RnBeads/reports/rnbSet_preprocessed", temp.dir = tempdir())




library("readxl")
SDRF <- read.csv("/Users/mariacastillo/Desktop/reports_ALL/differential_methylation_data/diffMethTable_S3.csv")
ensids <- SDRF$symbol
cols <- c('SYMBOL',"GENENAME")
gene_names=select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="SYMBOL")
colnames(gene_names) <- c('symbol',"GENENAME")
aa=merge(SDRF , gene_names, by='symbol')
write.csv(aa, "/Users/mariacastillo/Desktop/diffMethTable_S3.csv")

enrich.go=load("/Users/mariacastillo/Desktop/HEATSTROKE/Methylation/results_final/reports_ALL/differential_methylation_data/differential_rnbDiffMeth/enrich.go.RData")


rnb.options("filtering.greedycut"=TRUE)
rnb.options("assembly" = "hg38")
# Disable intersample variation plots (exploratory analysis)
rnb.options("exploratory.intersample"=TRUE)
# Reduce the subsampling number for estimating density plots
rnb.options("distribution.subsample"=10000)
# Disable regional methylation profiling (exploratory analysis)
rnb.options("exploratory.region.profiles"=NULL)
rnb.options("identifiers.column"="sampleID")
rnb.options("import.bed.style" = "bismarkCov")

# Disable chromosome coverage plots (QC, sequencing data only)
rnb.options("qc.boxplots"=TRUE)
rnb.options("qc.barplots"=TRUE)
rnb.options("qc.negative.boxplot"=TRUE)
rnb.options("qc.coverage.plots"=FALSE)
rnb.options("qc.coverage.histograms"=TRUE)
rnb.options("qc.coverage.violins"=TRUE)
rnb.options("qc.sample.batch.size"=500)

rnb.options("normalization" = TRUE)
rnb.options("normalization.plot.shifts"=TRUE)
rnb.options("filtering.greedycut"=FALSE)
rnb.options("filtering.greedycut.pvalue.threshold"=0.05)
rnb.options("filtering.low.coverage.masking" = FALSE)
rnb.options("filtering.high.coverage.outliers" = FALSE)
rnb.options("filtering.sex.chromosomes.removal"=TRUE)
rnb.options("filtering.deviation.threshold" = 0)
rnb.options("filtering.missing.value.quantile" = 0.5)
rnb.options("imputation.method" = "mean.samples")

rnb.options("exploratory.columns" = c("differentiation_level_1","Predicted Sex","Age"))
rnb.options("exploratory.clustering.heatmaps.pdf" = TRUE)
rnb.options("exploratory.principal.components"=8)
rnb.options("exploratory.correlation.pvalue.threshold"=0.01)
rnb.options("exploratory.correlation.permutations"=10000)
rnb.options("exploratory.correlation.qc"=TRUE)
rnb.options("exploratory.beta.distribution"=TRUE)
rnb.options("exploratory.intersample"=TRUE)
rnb.options("exploratory.deviation.plots"=TRUE)
rnb.options("exploratory.clustering.top.sites"=1000)

rnb.options("differential"=TRUE)
rnb.options("differential.comparison.columns"=c("differentiation_level_1","Predicted Sex","Age") )
rnb.options("differential.variability" = TRUE)
rnb.options("differential.variability.method" = "diffVar")
rnb.options("differential.enrichment.go"=TRUE)
rnb.options("differential.report.sites" = TRUE)
rnb.options("export.to.bed" = TRUE)
rnb.options("export.to.csv"=TRUE)
export.types = c("sites","tiling","genes","promoters","cpgislands")





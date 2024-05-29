library(methylKit)
library(ChIPseeker)

setwd("/Users/mariacastillo/Desktop/HEATSTROKE/CpG")
file_paths=read.delim("/Users/mariacastillo/Desktop/HEATSTROKE/CpG/files.txt")
samples=read.delim("/Users/mariacastillo/Desktop/HEATSTROKE/CpG_OT/filenames.txt")
infile <- system.file(list(file_paths), package = 'bsseq')
methylation <- read.bismark(file_paths, sample.ids = c(samples))
# Filter the CpG sites and normalize the data using the filterByCoverage and normalizeCoverage functions
methylation_filtered <- filterByCoverage(methylation, coverage = 10, lo.count = 3)
methylation_normalized <- normalizeCoverage(methylation_filtered)

# Perform differential methylation analysis using the calculateDiffMeth and getMethylDiff functions
diffmeth <- calculateDiffMeth(methylation_normalized, group1 = c("sample1"), group2 = c("sample2"))
results <- getMethylDiff(diffmeth, difference = 25, qvalue = 0.05)


peakAnno <- annotatePeak(results, tssRegion = c(-5000, 5000), TxDb = "org.Hs.eg.db", annoDb = "org.Hs.eg.db")

infile <- system.file(file_paths, package = 'methylKit')
myobj=methRead(read.bismark(file_paths, sample.id=list(samples),assembly="hg38"))



file.list=list(
  system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
  system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
  system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
  system.file("extdata", "control2.myCpG.txt", package = "methylKit") 
)         
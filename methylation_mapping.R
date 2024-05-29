#______________loading packages__________________#
library(TCGAbiolinks)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(RColorBrewer)
library(edgeR)



library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
man19 <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
grman19 <- makeGRangesFromDataFrame(data.frame(chr=man19$chr,
                                               start=man19$pos,
                                               end=man19$pos,
                                               cg.id=man19$Name,
                                               stringsAsFactors = F),
                                    keep.extra.columns = T)
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(path)
grman38 <- liftOver(grman19, ch)
grman38 <- unlist(grman38)
length(grman38) # 485344, = 485344/485512 = 99.97%, 168 CGs not re-mapped/dropped by liftOver

# subset, sort, and assemble hg38 manifest with hg19 position info included
grman19.lo <- grman19[grman19$cg.id %in% grman38$cg.id]
identical(grman19.lo$cg.id,grman38$cg.id)
identical(seqnames(grman19.lo),seqnames(grman38))

man38 <- data.frame(chr=seqnames(grman38),
                    pos=start(grman38),
                    Name=grman38$cg.id,
                    chr.hg19=seqnames(grman19.lo),
                    pos.hg19=start(grman19.lo),
                    stringsAsFactors=F)
rownames(man38) <- man38$Name

# subset, sort, assemble the gene and islands annotation from hg19
man19.lo <- man19[man19$Name %in% man38$Name,]
identical(man19.lo$Name,man38$Name)
man38$strand.hg19 <- man19.lo$strand
man38$ucsc.refgene.group.hg19 <- man19.lo$UCSC_RefGene_Group
man38$ucsc.refgene.name.hg19 <- man19.lo$UCSC_RefGene_Name
man38$cpg.island.name.hg19 <- man19.lo$Islands_Name
man38$cpg.island.region.hg19 <- man19.lo$Relation_to_Island



















#increase memory size
memory.limit(size = 28000) # do this on your own risk! 

#_________ DNA methylation Download____________#
# DNA methylation aligned to hg19
query_met <- GDCquery(project= "TCGA-BLCA", 
                      data.category = "DNA methylation", 
                      platform = "Illumina Human Methylation 450", 
                      legacy = TRUE)
GDCdownload(query_met)
#putting files togathers
data.met <- GDCprepare(query_met)
#saving the met object
saveRDS(object = data.met,
        file = "data.met.RDS",
        compress = FALSE)
# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.met.RDS")
# met matrix
met <- as.data.frame(SummarizedExperiment::assay(data.met))
# clinical data
clinical <- data.frame(data.met@colData)

#___________inspectiong methylation data_______________#

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## remove probes with NA
probe.na <- rowSums(is.na(met))

table(probe.na == 0)
#FALSE   TRUE 
#103553 382024 
# chose those has not NA values in rows
probe <- probe.na[probe.na == 0]
met <- met[row.names(met) %in% names(probe), ]

## remove probes that match to chromosome  X and Y 
keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
met <- met[keep, ]
rm(keep) # remove no further needed probes.

## remove SNPs overlapped probe
table(is.na(ann450k$Probe_rs))
# probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]

snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]

# filtre met
met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]

#remove no-further needed dataset
rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)

## Removing probes that have been demonstrated to map to multiple places in the genome.
# list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470

crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]

# filtre met
met <- met[ -which(row.names(met) %in% crs.reac), ]
bval <- met

## converting beta values to m_values
## m = log2(beta/1-beta)
mval <- t(apply(met, 1, function(x) log2(x/(1-x))))
#______________saving/loading_____________________#
# save data sets
#saveRDS(mval, file = "mval.RDS", compress = FALSE)
#saveRDS (bval, file = "bval.RDS", compress = FALSE)
#mval <- readRDS("mval.RDS")
#bval <- readRDS("bval.RDS")
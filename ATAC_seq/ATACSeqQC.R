library(BiocManager)
#BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
#                       "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene",
#                       "phastCons100way.UCSC.hg19"))

library(ATACseqQC)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GO.db)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
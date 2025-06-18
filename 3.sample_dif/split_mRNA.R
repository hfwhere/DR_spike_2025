library(tidyverse)
library(rtracklayer)
gtf <- as.data.frame(
  rtracklayer::import("genes.gtf")) %>%
  dplyr::select(transcript_id, transcript_biotype) %>%
  distinct() 

trans_counts <- read.table(file = "counts.matrix")
trans_exp <- read.table(file = "TMM.EXPR.matrix")

mrna_id = gtf[gtf$transcript_biotype == 'protein_coding','transcript_id']
lncrna_id = gtf[gtf$transcript_biotype == 'lncRNA','transcript_id']

write.table(trans_counts[mrna_id,], file = "mRNA.counts.matrix",
            sep = '\t',
            quote = F)

write.table(trans_exp[mrna_id,], file = "mRNA.TMM.EXPR.matrix",
            sep = '\t',
            quote = F)

write.table(trans_counts[lncrna_id,], file = "lncRNA.counts.matrix",
            sep = '\t',
            quote = F)

write.table(trans_exp[lncrna_id,], file = "lncRNA.TMM.EXPR.matrix",
            sep = '\t',
            quote = F)

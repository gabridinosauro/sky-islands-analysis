# Load required packages
# Set library path
.libPaths("/contrib/laurenpetrullo/dada2/lib/R/library")

library("dada2")
library("ggpubr")

setwd("/xdisk/laurenpetrullo/schiro/vsearch_GTDB")
#seqtab.nochim = readRDS("seqtab_nochim.rds")
#asvtab = t(seqtab.nochim)

### taxonomy
#seqs = colnames(asvtab)
#taxaGTDBfull <- assignTaxonomy(seqs, "/xdisk/laurenpetrullo/schiro/vsearch_GTDB/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz", multithread=TRUE)
#saveRDS(taxaGTDBfull,"taxaGTDBfull.RDS")


ref <- "/xdisk/laurenpetrullo/schiro/vsearch_GTDB/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz"
asvtab <- t(readRDS("seqtab_nochim.rds"))
seqs <- colnames(asvtab)
chunks <- split(seqs, ceiling(seq_along(seqs)/20000))
res <- lapply(chunks, function(ch) assignTaxonomy(ch, ref, multithread=TRUE))
taxa <- do.call(rbind, res)
saveRDS(taxa, "taxaGTDBfull.RDS")
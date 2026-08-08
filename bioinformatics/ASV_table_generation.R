library(Biostrings)
source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")

setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/all_data_run/run_dada2/")

#### Read in raw data
asvtab = readRDS("seqtab_nochim.rds")
taxtab = readRDS("taxaGTDBfull.RDS")

#### Build named ASV table + taxonomy table
asvtab_t = data.frame(t(asvtab))
sets = ASV_name_set(asvtab_t, taxtab)
ASVtable = sets$ASVtab
Taxtab   = as.data.frame(sets$TAXtab)

#### Save both tables
saveRDS(Taxtab, "Taxtab.RDS")
saveRDS(ASVtable, "ASVtable.RDS")

#### Write fasta of all ASVs (name -> sequence)
fasta_lines <- paste0(">", colnames(ASVtable), "\n", colnames(asvtab_t))
writeLines(fasta_lines, "sorted_ASVs.fasta")

#### Filter ASV table (min 50 total reads, min 10 reads in at least one sample)
ASVtable_filt <- ASVtable[, colSums(ASVtable) >= 50 & apply(ASVtable, 2, max) >= 10]

#### Subset sequences to filtered ASVs, write final fasta
seqs <- readDNAStringSet("sorted_ASVs.fasta")
seqs1 <- seqs[names(seqs) %in% colnames(ASVtable_filt)]
writeXStringSet(seqs1, filepath = "filtered_ASVs.fasta")
#### 01 - Mantel test computation (HPC-heavy step) ##########################
# Loops over every node subtree, computes Spearman & Pearson Mantel tests
# between ASV phylogenetic distance and host phylogenetic distance.
# This is the computationally expensive part - saves a checkpoint RDS
# ("result_final.RDS") that script 02 picks up from.
###############################################################################

library(vegan)
library(ape)
library(dplyr)

setwd("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")

## ---- Load inputs ----------------------------------------------------------
tree_list    <- list.files("/Users/gabri/Desktop/mouse/SILVA_tree/subtrees", full.names = TRUE)
tree_host    <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/expanded_tree.RDS")
ASVtab_rar   <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")
Taxtab       <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/Taxtab.RDS")
metadata_rar <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/metadata.RDS")

host_phy_mat <- cophenetic.phylo(tree_host)

## ---- Chunked Mantel loop ---------------------------------------------------
results_final <- data.frame(
  node = NULL, mantel_stat_spear = NULL, mantel_p_spear = NULL,
  mantel_stat_pear = NULL, mantel_p_pear = NULL
)

chunk_size <- 100
chunks <- split(seq_along(tree_list), ceiling(seq_along(tree_list) / chunk_size))

for (chunk_idx in seq_along(chunks)) {
  cat("Processing chunk", chunk_idx, "of", length(chunks), "\n")
  
  chunk_results <- data.frame(
    node = NULL, mantel_stat_spear = NULL, mantel_p_spear = NULL,
    mantel_stat_pear = NULL, mantel_p_pear = NULL
  )
  
  for (i in chunks[[chunk_idx]]) {
    tree_an   <- read.tree(tree_list[i])
    node_name <- tools::file_path_sans_ext(basename(tree_list[i]))
    
    tips <- tree_an$tip.label
    ASVtab_filt <- ASVtab_rar[, colnames(ASVtab_rar) %in% tips]
    ASVtab_filt <- ASVtab_filt[rowSums(ASVtab_filt) > 0, ]
    result <- data.frame(
      sample = rownames(ASVtab_filt),
      most_abundant_ASV = colnames(ASVtab_filt)[apply(ASVtab_filt, 1, which.max)],
      abundance = apply(ASVtab_filt, 1, max)
    )
    tips_to_keep <- tips[tips %in% result$most_abundant_ASV]
    tree_an_1 <- keep.tip(tree_an, tips_to_keep)
    ASV_phy_mat <- cophenetic.phylo(tree_an_1)
    
    shared_samples <- intersect(result$sample, rownames(host_phy_mat))
    result_filt <- result[result$sample %in% shared_samples, ]
    asv_names <- result_filt$most_abundant_ASV
    names(asv_names) <- result_filt$sample
    asv_dist_mat <- matrix(NA, nrow = nrow(result_filt), ncol = nrow(result_filt),
                           dimnames = list(result_filt$sample, result_filt$sample))
    for (si in result_filt$sample) {
      for (sj in result_filt$sample) {
        asv_dist_mat[si, sj] <- ASV_phy_mat[asv_names[si], asv_names[sj]]
      }
    }
    host_dist_mat <- host_phy_mat[result_filt$sample, result_filt$sample]
    mt_spear <- mantel(as.dist(asv_dist_mat), as.dist(host_dist_mat), method = "spearman",
                       permutations = 999, parallel = getOption("mc.cores"))
    mt_pear  <- mantel(as.dist(asv_dist_mat), as.dist(host_dist_mat), method = "pearson",
                       permutations = 999, parallel = getOption("mc.cores"))
    
    chunk_results <- rbind(chunk_results, data.frame(
      node = node_name,
      mantel_stat_spear = mt_spear$statistic, mantel_p_spear = mt_spear$signif,
      mantel_stat_pear  = mt_pear$statistic,  mantel_p_pear  = mt_pear$signif
    ))
  }
  
  # Checkpoint each chunk to disk (lets you resume/parallelize per-chunk on SLURM)
  saveRDS(chunk_results, file = paste0("mantel_chunk_", chunk_idx, ".rds"))
  cat("  Saved chunk", chunk_idx, "\n")
  
  results_final <- rbind(results_final, chunk_results)
}

## ---- Final checkpoint for script 02 ---------------------------------------
saveRDS(results_final, "result_final.RDS")
cat("Done. Saved result_final.RDS with", nrow(results_final), "nodes.\n")
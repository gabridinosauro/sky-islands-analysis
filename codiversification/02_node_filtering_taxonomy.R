#### 02 - Node filtering, redundancy removal & taxonomy annotation ##########
# Picks up from result_final.RDS (script 01).
# - Aligns Mantel results to node list
# - Removes redundant nodes (>=75% tip overlap, keeping the higher Mantel r)
# - Annotates each retained node with its lowest common taxonomy
# - FDR-corrects p-values, defines the significant/strong codiversifying set
# - Identifies dominant vs non-dominant ASVs per node (used for tree pruning
#   and plotting in script 03)
###############################################################################

library(dplyr)
library(ape)
library(pbapply)

setwd("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")

## ---- Load checkpoint from script 01 + node metadata ------------------------
results_final <- readRDS("result_final.RDS")
nodes_list    <- readRDS("/Users/gabri/Desktop/mouse/SILVA_tree/treeees/node_list_filtered.RDS")
Taxtab        <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/Taxtab.RDS")
ASVtab_rar    <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")

results_final$node_name <- gsub("subtree_node", "", results_final$node)
results_final <- results_final[results_final$node_name %in% nodes_list$node_id, ]

idx <- match(nodes_list$node_id, results_final$node_name)
results_final <- results_final[idx[!is.na(idx)], ]

idx <- match(nodes_list$node_id, results_final$node_name)
results_final_aligned <- results_final[idx, ]
unique(nodes_list$node_id == results_final_aligned$node_name) ## sanity check - should be TRUE

results_tot <- cbind(nodes_list, results_final_aligned)

## ---- Remove redundant nodes (keep highest-r node among overlapping sets) --
filter_redundant_nodes <- function(df, overlap_threshold = 0.75) {
  df <- df %>% arrange(desc(mantel_stat_spear))
  
  kept <- c()
  removed <- c()
  
  for (i in seq_len(nrow(df))) {
    focal_id <- df$node_id[i]
    if (focal_id %in% removed) next
    
    kept <- c(kept, focal_id)
    tips_focal <- unlist(df$descendants[i])
    
    for (j in seq_len(nrow(df))) {
      candidate_id <- df$node_id[j]
      if (candidate_id == focal_id || candidate_id %in% removed || candidate_id %in% kept) next
      
      tips_candidate <- unlist(df$descendants[j])
      n_shared <- length(intersect(tips_focal, tips_candidate))
      overlap <- n_shared / min(length(tips_focal), length(tips_candidate))
      
      if (overlap >= overlap_threshold) {
        removed <- c(candidate_id, removed)
      }
    }
  }
  df %>% filter(node_id %in% kept)
}

results_filtered <- filter_redundant_nodes(results_tot, overlap_threshold = 0.75)
results_filtered$ASV <- lapply(results_filtered$descendants, names)

## ---- Lowest common taxonomy per retained node ------------------------------
get_lowest_common_tax <- function(asvs, taxtab) {
  tax_sub <- taxtab[rownames(taxtab) %in% asvs, , drop = FALSE]
  if (nrow(tax_sub) == 0) return(NA_character_)
  
  levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  for (lvl in levels) {
    vals <- tax_sub[[lvl]]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0 && length(unique(vals)) == 1) {
      return(paste0(lvl, ": ", unique(vals)))
    }
  }
  return("No consensus")
}

results_filtered$lowest_common_tax <- sapply(
  results_filtered$ASV, get_lowest_common_tax, taxtab = Taxtab
)

results_export <- data.frame(
  lapply(results_filtered, function(x) if (is.list(x)) sapply(x, toString) else x),
  check.names = FALSE
)
#write.csv(results_export, "results_filtered.csv", row.names = FALSE)

## ---- FDR correction & significant/strong codiversifying set ---------------
results_filtered$p_cor <- p.adjust(results_filtered$mantel_p_spear, method = "fdr")

Mantel_res_sig      <- results_filtered[results_filtered$p_cor < 0.05, ]
Mantel_res_sig_filt <- Mantel_res_sig[Mantel_res_sig$mantel_stat_spear > 0.75, ]

all_asvs <- unique(unlist(lapply(Mantel_res_sig_filt$descendants, names)))
tax_sig  <- Taxtab[all_asvs, ]
print(sort(table(tax_sig$Family)))

## ---- Dominant vs non-dominant ASV per node (needed for tree pruning) ------
results <- pbsapply(nodes_list$descendants, function(desc) {
  asvs <- names(desc)[names(desc) %in% colnames(ASVtab_rar)]
  if (length(asvs) == 0) return(list(dominant = character(0), non_dominant = character(0)))
  
  sub <- ASVtab_rar[, asvs, drop = FALSE]
  present_rows <- rowSums(sub) > 0
  sub <- sub[present_rows, , drop = FALSE]
  if (nrow(sub) == 0) return(list(dominant = character(0), non_dominant = character(0)))
  
  dominant_asvs <- unique(colnames(sub)[apply(sub, 1, which.max)])
  non_dominant_asvs <- asvs[!asvs %in% dominant_asvs]
  list(dominant = dominant_asvs, non_dominant = non_dominant_asvs)
}, simplify = FALSE)

nodes_list$dominant_asvs     <- lapply(results, `[[`, "dominant")
nodes_list$non_dominant_asvs <- lapply(results, `[[`, "non_dominant")

## ---- Mean relative abundance of dominant ASVs per node ---------------------
ASVtab_rel <- ASVtab_rar / rowSums(ASVtab_rar)

nodes_list$dominant_mean_relabund <- sapply(nodes_list$dominant_asvs, function(asvs_str) {
  asvs <- trimws(unlist(strsplit(as.character(asvs_str), ",")))
  asvs <- asvs[asvs %in% colnames(ASVtab_rel)]
  if (length(asvs) == 0) return(NA)
  vals <- unlist(ASVtab_rel[, asvs, drop = FALSE])
  mean(vals[vals > 0])
})

Mantel_res_sig_filt <- Mantel_res_sig_filt %>%
  left_join(nodes_list %>% select(node_id, dominant_mean_relabund), by = "node_id")

## FDR-corrected version of the *unfiltered* node set (used again in script 04)
results_final$padj <- p.adjust(results_final$mantel_p_spear, method = "fdr")

## ---- Checkpoints for scripts 03 & 04 ---------------------------------------
saveRDS(nodes_list,          "nodes_list_1.RDS")
saveRDS(results_filtered,    "results_filtered_new.RDS")
saveRDS(Mantel_res_sig_filt, "Mantel_res_sig_filt.RDS")
saveRDS(results_final,       "results_final_padj.RDS")



  #### 05 - Mantel tests split by mountain range ################################
# Picks up Mantel_res_sig_filt.RDS from script 02.
# For every strongly codiversifying node (r > 0.75, padj < 0.05, from the
# combined dataset), re-runs the host-microbe Mantel test separately within
# each mountain range (Pinaleño vs Santa Catalina) to check whether the
# codiversification signal holds within a mountain or is mostly driven by
# between-mountain host divergence.
###############################################################################

library(vegan)
library(ape)
library(ggplot2)
library(patchwork)

setwd("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")

## ---- Load checkpoint (script 02) + raw inputs ------------------------------
Mantel_res_sig_filt <- readRDS("Mantel_res_sig_filt.RDS")
tree_list    <- list.files("/Users/gabri/Desktop/mouse/SILVA_tree/subtrees", full.names = TRUE)
ASVtab_rar   <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")
metadata_rar <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/metadata.RDS")
tree_host    <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/expanded_tree.RDS")

Mantel_res_sig_filt <- na.omit(Mantel_res_sig_filt)

## ---- Per-node, per-mountain Mantel tests -----------------------------------
mantel_results_df <- data.frame(
  node_id = character(nrow(Mantel_res_sig_filt)),
  pin_r   = numeric(nrow(Mantel_res_sig_filt)),
  pin_p   = numeric(nrow(Mantel_res_sig_filt)),
  sc_r    = numeric(nrow(Mantel_res_sig_filt)),
  sc_p    = numeric(nrow(Mantel_res_sig_filt))
)

for (i in seq_len(nrow(Mantel_res_sig_filt))) {
  name         <- Mantel_res_sig_filt$node_id[i]
  tree_file    <- tree_list[grep(name, tree_list)]
  tree_to_plot <- read.tree(tree_file)
  tips <- tree_to_plot$tip.label

  ASVtab_filt <- ASVtab_rar[, colnames(ASVtab_rar) %in% tips]
  ASVtab_filt <- ASVtab_filt[rowSums(ASVtab_filt) > 0, ]

  result <- data.frame(
    sample            = rownames(ASVtab_filt),
    most_abundant_ASV = colnames(ASVtab_filt)[apply(ASVtab_filt, 1, which.max)],
    abundance         = apply(ASVtab_filt, 1, max)
  )
  result$mountain <- metadata_rar$Mountain[match(result$sample, metadata_rar$sample_id)]

  ### Pinaleño
  result_pin          <- result[result$mountain == "Pinaleno", ]
  tips_to_keep         <- tips[tips %in% result_pin$most_abundant_ASV]
  tree_pin             <- keep.tip(tree_to_plot, tips_to_keep)
  final_tree_host_pin  <- keep.tip(tree_host, result_pin$sample)

  dist_micro_pin <- cophenetic(tree_pin)[result_pin$most_abundant_ASV, result_pin$most_abundant_ASV]
  dist_host_pin  <- cophenetic(final_tree_host_pin)[result_pin$sample, result_pin$sample]

  mantel_pin <- vegan::mantel(as.dist(dist_micro_pin), as.dist(dist_host_pin),
                        method = "spearman", permutations = 999)

  ### Santa Catalina
  result_sc          <- result[result$mountain == "Santa Catalina", ]
  tips_to_keep        <- tips[tips %in% result_sc$most_abundant_ASV]
  tree_sc             <- keep.tip(tree_to_plot, tips_to_keep)
  final_tree_host_sc  <- keep.tip(tree_host, result_sc$sample)

  dist_micro_sc <- cophenetic(tree_sc)[result_sc$most_abundant_ASV, result_sc$most_abundant_ASV]
  dist_host_sc  <- cophenetic(final_tree_host_sc)[result_sc$sample, result_sc$sample]

  mantel_sc <- vegan::mantel(as.dist(dist_micro_sc), as.dist(dist_host_sc),
                       method = "spearman", permutations = 999)

  mantel_results_df[i, ] <- list(name,
                                  mantel_pin$statistic, mantel_pin$signif,
                                  mantel_sc$statistic,  mantel_sc$signif)
  cat("Node", i, "of", nrow(Mantel_res_sig_filt), "done\n")
}

## ---- Plot: distribution of within-mountain Mantel r, colored by significance
mantel_results_df$pin_sig <- factor(
  ifelse(mantel_results_df$pin_p < 0.05, "significant", "non-significant"),
  levels = c("non-significant", "significant")
)
mantel_results_df$sc_sig <- factor(
  ifelse(mantel_results_df$sc_p < 0.05, "significant", "non-significant"),
  levels = c("non-significant", "significant")
)

sig_colors <- c("non-significant" = "grey70", "significant" = "#2c7fb8")

max_count <- max(
  hist(mantel_results_df$pin_r, breaks = seq(0, 1, length.out = 51), plot = FALSE)$counts,
  hist(mantel_results_df$sc_r,  breaks = seq(0, 1, length.out = 51), plot = FALSE)$counts
)

p1 <- ggplot(mantel_results_df, aes(x = pin_r, fill = pin_sig)) +
  geom_histogram(bins = 50, color = "white") +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "black") +
  scale_fill_manual(values = sig_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, max_count)) +
  labs(title = "Pinaleño", x = "Mantel rho", y = "Count", fill = NULL) +
  theme_classic()

p2 <- ggplot(mantel_results_df, aes(x = sc_r, fill = sc_sig)) +
  geom_histogram(bins = 50, color = "white") +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "black") +
  scale_fill_manual(values = sig_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, max_count)) +
  labs(title = "Santa Catalina", x = "Mantel rho", y = "Count", fill = NULL) +
  theme_classic()

(p1 / p2) + plot_layout(guides = "collect", axes = "collect")

## ---- Checkpoint for downstream use -----------------------------------------
saveRDS(mantel_results_df, "mantel_results_by_mountain.RDS")
cat("Done. Saved mantel_results_by_mountain.RDS (", nrow(mantel_results_df), "nodes).\n")

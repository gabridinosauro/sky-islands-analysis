#### 03 - Tree visualization ##################################################
# Picks up from checkpoints saved in script 02.
# Part A: single-node cophylo plot (host tree vs one ASV subtree, e.g. node 70768)
# Part B: whole-tree circular ggtree colored by Mantel's r, with a phylum ring
###############################################################################

library(dplyr)
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(tidytree)
library(phytools)
library(RColorBrewer)
library(viridisLite)
library(pbapply)

setwd("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")

## ---- Load checkpoints & raw data -------------------------------------------
nodes_list          <- readRDS("nodes_list_1.RDS")
Mantel_res_sig_filt  <- readRDS("Mantel_res_sig_filt.RDS")
results_filtered     <- readRDS("results_filtered_new.RDS")

tree_list    <- list.files("/Users/gabri/Desktop/mouse/SILVA_tree/subtrees", full.names = TRUE)
tree_host    <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/expanded_tree.RDS")
ASVtab_rar   <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")
Taxtab       <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/Taxtab.RDS")
metadata_rar <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/metadata.RDS")
tree_ASV     <- readRDS("/Users/gabri/Desktop/mouse/SILVA_tree/tt_no_archaea.RDS")

## =============================================================================
## PART A - single-node cophylo plot (example node "70768")
## Swap the grep pattern to plot a different node id.
## =============================================================================
tree_file    <- tree_list[grep("70768", tree_list)]
tree_to_plot <- read.tree(tree_file)

tips <- tree_to_plot$tip.label
ASVtab_filt <- ASVtab_rar[, colnames(ASVtab_rar) %in% tips]
ASVtab_filt <- ASVtab_filt[rowSums(ASVtab_filt) > 0, ]
result <- data.frame(
  sample = rownames(ASVtab_filt),
  most_abundant_ASV = colnames(ASVtab_filt)[apply(ASVtab_filt, 1, which.max)],
  abundance = apply(ASVtab_filt, 1, max)
)
tips_to_keep <- tips[tips %in% result$most_abundant_ASV]
tree_to_plot <- keep.tip(tree_to_plot, tips_to_keep)

ass <- data.frame(V1 = result$sample, V2 = result$most_abundant_ASV)
ass1 <- ass
ass1$species <- metadata_rar$species[match(ass1$V1, metadata_rar$sample_id)]

species_colors <- c(
  "Peromyscus melanotis" = "#fb6a4a", "Peromyscus boylii" = "#de2d26",
  "Peromyscus eremicus" = "#a50f15", "Neotoma albigula" = "#fd8d3c",
  "Neotoma mexicana" = "#e6550d", "Microtus longicaudus" = "#fec44f",
  "Reithrodontomys megalotis" = "#fee391", "Sigmodon ochrognathus" = "#d4a017",
  "Onychomys torridus" = "#c51b8a", "Chaetodipus penicillatus" = "#6baed6",
  "Chaetodipus baileyi" = "#3182bd", "Chaetodipus intermedius" = "#9ecae1",
  "Dipodomys merriami" = "#08519c", "Neotamias dorsalis" = "#41ab5d",
  "Ammospermophilus harrisii" = "#006d2c"
)

ass1$colors <- species_colors[match(ass1$species, names(species_colors))]
final_tree_host <- keep.tip(tree_host, ass1$V1)

mountain_colors <- c("Pinaleno" = "red", "Santa Catalina" = "blue")
cophylo_obj <- cophylo(final_tree_host, tree_to_plot, assoc = ass, rotate = TRUE)

plot(cophylo_obj,
     type = "phylogram", link.type = "curved", link.lty = "solid",
     link.lwd = 1, link.col = ass1$colors, ftype = "off", fsize = 0.3, cex = 0.5)

left_tips    <- cophylo_obj$trees[[1]]$tip.label
tip_mountain <- metadata_rar$Mountain[match(left_tips, metadata_rar$sample_id)]
left_cols    <- unname(mountain_colors[tip_mountain])
left_cols[is.na(left_cols)] <- "black"

tiplabels.cophylo(which = "left", pch = 21, bg = left_cols, col = left_cols, cex = 1)

# Mantel test restricted to this single node's tree (sanity check against the loop result)
final_dist_host <- cophenetic.phylo(final_tree_host)
ASV_phy_mat <- cophenetic.phylo(keep.tip(tree_to_plot, tips_to_keep))
shared_samples <- intersect(result$sample, rownames(final_dist_host))
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
asv_dist_mat <- asv_dist_mat[rownames(final_dist_host), colnames(final_dist_host)]

mt_spear <- vegan::mantel(as.dist(asv_dist_mat), as.dist(final_dist_host),
                           method = "spearman", permutations = 999, parallel = getOption("mc.cores"))
mt_spear$statistic
mt_spear$signif

plot_df <- data.frame(
  Host_distance = as.vector(as.dist(final_dist_host)),
  ASV_distance  = as.vector(as.dist(asv_dist_mat))
)
ggplot(plot_df, aes(x = Host_distance, y = ASV_distance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_classic() +
  labs(x = "Host phylogenetic distance", y = "ASV phylogenetic distance")

## =============================================================================
## PART B - whole-tree circular ggtree colored by Mantel's r + phylum ring
## =============================================================================
all_dominant <- unique(unlist(nodes_list$dominant_asvs, use.names = FALSE))

Ntip_full  <- Ntip(tree_ASV)
Nnode_full <- tree_ASV$Nnode
node_label_vals <- rep(NA_real_, Nnode_full)

for (i in seq_len(nrow(Mantel_res_sig_filt))) {
  full_node_id <- Mantel_res_sig_filt$node_id[i]
  idx <- full_node_id - Ntip_full
  if (!is.na(idx) && idx >= 1 && idx <= Nnode_full) {
    node_label_vals[idx] <- Mantel_res_sig_filt$mantel_stat_spear[i]
  }
}
tree_ASV$node.label <- as.character(node_label_vals)

pruned_tree <- keep.tip(tree_ASV, all_dominant)
saveRDS(pruned_tree, "pruned_tree.RDS")

node_vals <- as.numeric(pruned_tree$node.label)
tree_tbl <- as_tibble(pruned_tree)

get_tip_value <- function(tree, tip_node, node_vals) {
  ancestors <- c()
  current <- tip_node
  while (!is.na(current)) {
    ancestors <- c(ancestors, current)
    current <- tree$edge[tree$edge[, 2] == current, 1]
    if (length(current) == 0) break
  }
  internal_ancestors <- ancestors[ancestors > Ntip(tree)]
  idx <- internal_ancestors - Ntip(tree)
  valid <- idx[!is.na(node_vals[idx])]
  if (length(valid) == 0) return(NA)
  return(node_vals[valid[1]])
}

tip_nodes  <- 1:Ntip(pruned_tree)
tip_values <- pbsapply(tip_nodes, function(tip) get_tip_value(pruned_tree, tip, node_vals))
saveRDS(tip_values, "tip_values.RDS")

tree_tbl$mantel_final <- NA_real_
tree_tbl$mantel_final[tree_tbl$node %in% tip_nodes] <- tip_values

p <- ggtree(pruned_tree, layout = "circular", branch.length = "none",
            aes(color = mantel_final), size = 0.4) %<+% tree_tbl
p + scale_color_viridis_c(option = "magma", na.value = "grey90", name = "Mantel's r")

## Subsample to 10,000 tips for a cleaner circular cladogram
set.seed(123)
n_keep <- min(10000, length(pruned_tree$tip.label))
random_tips <- sample(pruned_tree$tip.label, n_keep)
subset_tree <- keep.tip(pruned_tree, random_tips)
subset_tbl  <- as_tibble(subset_tree)

if (!("isTip" %in% names(tree_tbl))) {
  tree_tbl <- tree_tbl %>% mutate(isTip = node <= Ntip(pruned_tree))
}
tip_values_df <- tree_tbl %>%
  filter(.data$isTip == TRUE) %>%
  select(label, mantel_final) %>%
  distinct()

subset_tbl <- subset_tbl %>% left_join(tip_values_df, by = "label")

p <- ggtree(subset_tree, layout = "circular", branch.length = "none", size = 0.02) %<+% subset_tbl +
  geom_tree(aes(color = mantel_final), size = 0.8, show.legend = TRUE) +
  scale_color_gradient(low = "#FEE5D9", high = "#A50F15", na.value = "grey95", name = "Mantel's r")
p

## ---- Add a phylum ring around the subsampled tree --------------------------
tax_df <- Taxtab %>% tibble::rownames_to_column("label") %>% select(label, Phylum)

ring_df <- data.frame(label = subset_tree$tip.label) %>%
  left_join(tax_df, by = "label") %>%
  rename(group = Phylum)

phylum_levels <- sort(unique(ring_df$group))
phylum_colors <- setNames(viridis(n = length(phylum_levels), option = "turbo"), phylum_levels)

p <- ggtree(subset_tree, layout = "circular", branch.length = "none",
            size = 0.02, lineend = "butt") %<+% subset_tbl +
  geom_tree(color = "grey90", size = 0.4) +
  geom_tree(aes(color = mantel_final, subset = !is.na(mantel_final)), size = 1.2) +
  scale_color_gradient(low = "gold", high = "#A50F15", na.value = NA, name = "Mantel's r") +
  geom_fruit(data = ring_df, geom = geom_tile,
             mapping = aes(y = label, fill = group), width = 1, offset = 0.1) +
  scale_fill_manual(values = phylum_colors, name = "Phylum") +
  geom_tippoint(aes(subset = (isTip & !is.na(mantel_final))), color = "red", size = 1.5)
p

## Checkpoint for script 04 (family-level summary needs tip-level taxonomy + mantel status)
saveRDS(tip_values_df, "tip_values_df.RDS")
cat("Done with tree visualization.\n")

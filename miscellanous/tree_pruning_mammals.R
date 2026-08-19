library(ape)
library(treeio)

setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/dada2/")
metadata <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/metadata_clean.RDS")
mam_full_tree <- treeio::read.beast("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")
phyphy = mam_full_tree@data

# Fix taxonomy and convert to tree-style labels
species_vec <- unique(metadata$species)
species_vec <- gsub("^Neotamias dorsalis$", "Tamias dorsalis", species_vec)
species_underscore <- gsub(" ", "_", species_vec)
species_vec
# Extract bare phylo FIRST, then prune with ape
phy_full <- as.phylo(mam_full_tree)
tip_labels <- phy_full$tip.label
phy_full$node.label <- mam_full_tree@data$posterior[order(mam_full_tree@data$node)]
pattern <- paste0("^_?(", paste(species_underscore, collapse = "|"), ")")
keep_tips <- grep(pattern, tip_labels, value = TRUE)


# Use ape::drop.tip - guaranteed to preserve branch lengths
phy <- ape::drop.tip(phy_full, setdiff(tip_labels, keep_tips))

# Clean tip labels
phy$tip.label <- sub("^_?([^_]+_[^_]+).*$", "\\1", phy$tip.label)

# Verify branch lengths survived
summary(phy$edge.length)

plot(phy, use.edge.length = TRUE)
add.scale.bar()
nodelabels(round(as.numeric(phy$node.label), 2), cex = 0.6)

saveRDS(phy, "Nates_tree.RDS")

saveRDS(phy, "Nates_tree.RDS")
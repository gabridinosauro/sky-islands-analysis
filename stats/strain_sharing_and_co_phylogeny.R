### Shared ASVs
### Packages -------------------------------------------------------------
library(vegan)
library(ecodist)
library(ape)
library(phytools)
library(reshape2)
library(stringr)
library(pvclust)
library(ggplot2)
library(ggpubr)
library(viridis)
library(cultevo)
source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")

setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/")
bac.dist.sor = readRDS("bac.dist_sor.RDS")
metadata = readRDS("metadata_clean_1.RDS")
ASVtable = readRDS("ASVtable.RDS")
Taxatable= readRDS("Taxtab.RDS")


### Check if samples in both tables are in the same order
unique(rownames(ASVtable) == metadata$sample_id) # OK
ASVtable = ASVtable[metadata$species != "Sorex monticolus", ] ### remove those samples from the ASvtable
metadata = metadata[metadata$species != "Sorex monticolus",]
unique(rownames(ASVtable) == metadata$sample_id) # OK

### read ASVs in--------------
mapping_file = readRDS("/Users/gabri/Desktop/mouse/SILVA_tree/nodes_df.RDS")
mapping_file = mapping_file[mapping_file$is_transition == TRUE,]
mapping_file$descendants[1]
all_asvs <- names(unlist(mapping_file$descendants))
ASVtable1 = ASVtable[,colnames(ASVtable) %in% all_asvs]
sort(rowSums(ASVtable1))
sort(rowSums(ASVtable))
### both 4 under 10000 --- ok
ASVtable = ASVtable1

#### remove samples with low count
sort(rowSums(ASVtable))
metadata$read_number = rowSums(ASVtable)
metadata = metadata[metadata$read_number > 10000, ] ### remove those samples from the metadata
ASVtable = ASVtable[rowSums(ASVtable) >10000, ] ### remove those samples from the ASvtable
ndvi_mod <- lm(mean_ndvi_2022_2025 ~ elevation, data = metadata)
metadata$ndvi_resid <- residuals(ndvi_mod)

###taxtab

taxtab_fin = Taxatable[rownames(Taxatable) %in% colnames(ASVtable), ]
### What is on average the share composition of members of the same species
names(bac.dist.sor)
ordered_tips = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/sky_islands_analysis/ordered_tips.RDS")



dm <- as.matrix(bac.dist.sor)
sp <- metadata$species[match(rownames(dm), metadata$sample_id)]
spp <- unique(sp)
m <- matrix(NA, nrow=length(spp), ncol=length(spp), dimnames=list(spp,spp))
for(i in seq_along(spp)) for(j in seq_along(spp))
  m[i,j] <- mean(dm[sp==spp[i], sp==spp[j]], na.rm=TRUE)
m
# Spend most colors on blue-green transition
cols <- colorRampPalette(c("darkblue","green","yellowgreen",
                           "greenyellow","yellow","orange","darkorange","red"),
                         bias = 1.2)(200)  # bias < 1 stretches low end
ord <- ordered_tips#hclust(as.dist(m))$order
m_ord <- m[ord, ord]

# Convert distance to similarity
m2 <- 1 - m_ord

# Compress everything below similarity 0.6

library(corrplot)
corrplot(m2,
         method = "color",
         type = "lower",
         is.corr = FALSE,
         diag = FALSE,
         col = cols,
         #col.lim = range(m2, na.rm = TRUE),
         tl.col = "black",
         tl.cex = 0.8)


### Collapse ASV table by species ----------------------------------------
# One row per species; no Mountain split
group_var <- metadata$species
asv_by_taxon <- rowsum(ASVtable, group = group_var)

# Check dimensions
nrow(asv_by_taxon)
ncol(asv_by_taxon)

### Match taxonomy table to ASV columns ---------------------------------
Taxtab1 <- Taxatable[rownames(Taxatable) %in% colnames(asv_by_taxon), ]

ASVtable_tax <- generate.tax.summary.modified(asv_by_taxon, Taxtab1)

### Phylum and family summaries -----------------------------------------
phylum_tab <- ASVtable_tax$tax3
phylum_values <- rowSums(phylum_tab)
phylum_values <- phylum_values / sum(phylum_values)
sort(phylum_values)

fam_tab <- ASVtable_tax$tax6
head(fam_tab)

family_values <- rowSums(fam_tab)
family_values <- family_values / sum(family_values)
sort(family_values)

### Family barplot --------------------------------------------------------
top10 <- names(sort(rowSums(fam_tab), decreasing = TRUE)[1:10])

fam_tab$Family <- rownames(fam_tab)
ASVtab_rel_melt <- reshape2::melt(fam_tab)
ASVtab_rel_melt$Family <- ifelse(ASVtab_rel_melt$Family %in% top10, ASVtab_rel_melt$Family, "others")

mypalette2 <- c(
  "#40004b","#ffffbf","#762a83","#de77ae","#f46d43","#f7f7f7","#d9f0d3","#a6dba0",
  "#a6cee3","#1f78b4","#b2df8a","#33a02c","#4d9221","#2166ac","#5aae61","#1b7837",
  "#d73027","#00441b","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#b2182b",
  "#de77ae","#bc80bd","#f5f5f5","#c7eae5","#80cdc1","#35978f","#003c30","#8e0152",
  "#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419",
  "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6",
  "#6a3d9a","#ffff99"
)

set.seed(3)
unique_cols <- unique(mypalette2)
mypalette20 <- sample(unique_cols, 12)

# Species order for plotting
species_order <- sort(unique(group_var))
ASVtab_rel_melt$variable <- factor(ASVtab_rel_melt$variable, levels = species_order)

ASVtab_rel_melt$variable_lab <- ASVtab_rel_melt$variable |>
  str_replace_all("_", " ")

ggplot(ASVtab_rel_melt, aes(fill = Family, y = value, x = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = mypalette20) +
  scale_x_discrete(labels = function(x) {
    parse(text = vapply(x, function(lbl) {
      paste0("italic('", lbl, "')")
    }, character(1)))
  }) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top"
  ) +
  ylab("Relative abundance") +
  xlab("Species") +
  guides(fill = guide_legend(nrow = 4))

### ---- SPLICED FROM SCRIPT 2: build tree_exp -----------------------------
### (needed below for tree reordering, phylogenetic distance, and cophylo)
tree_exp = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/dada2/Nates_tree.RDS")
is.ultrametric(tree_exp)
is.rooted(tree_exp)
beast_nwk <- write.tree(tree_exp)
tree_exp <- read.tree(text = beast_nwk)
is.rooted(tree_exp)

tree_exp = drop.tip(tree_exp, "Sorex_monticolus")
tree_exp$tip.label[tree_exp$tip.label == "Tamias_dorsalis"] = "Neotamias_dorsalis"
### Underscores -> spaces so tip labels match metadata$species / species_order
tree_exp$tip.label <- gsub("_", " ", tree_exp$tip.label)
### -------------------------------------------------------------------------

### Tree reordering ------------------------------------------------------
species_order_plot <- rev(species_order)

# If your tree labels have underscores, clean them here
tree_exp$tip.label <- gsub("_", " ", tree_exp$tip.label)

tree_exp2 <- rotateConstr(tree_exp, species_order_plot)
plot(tree_exp2)

### Rarefaction and average Bray distance --------------------------------
set.seed(42)
nrep <- 100
depth <- min(10000, min(rowSums(asv_by_taxon)))

dist_list <- vector("list", nrep)

for (i in seq_len(nrep)) {
  rar <- rrarefy(asv_by_taxon, sample = depth)
  dist_list[[i]] <- as.matrix(vegdist(rar, method = "bray", binary = TRUE))
}

avg_dist_mat <- Reduce(`+`, dist_list) / length(dist_list)

### ---- SPLICED FROM SCRIPT 2: build tree_avg from avg_dist_mat -----------
hc_avg <- hclust(as.dist(avg_dist_mat), method = "average")
tree_avg <- as.phylo(hc_avg)
plot(tree_avg)
### -------------------------------------------------------------------------

### Phylogenetic distance matrix -----------------------------------------
phy_mat1 <- cophenetic.phylo(tree_exp)

# Keep only common species
common <- intersect(rownames(phy_mat1), rownames(avg_dist_mat))
phy_mat1 <- phy_mat1[common, common]
avg_dist_mat <- avg_dist_mat[common, common]

# Convert to distance objects
phy_dist <- as.dist(phy_mat1)
micro_dist <- as.dist(avg_dist_mat)

### Mantel test ----------------------------------------------------------
mantel_result <- vegan::mantel(
  phy_dist,
  micro_dist,
  method = "spearman",
  permutations = 9999
)

mantel_result$statistic
mantel_result$signif

### Scatter plot ---------------------------------------------------------
df <- data.frame(
  phy = as.vector(phy_dist),
  micro = as.vector(micro_dist)
)

ggplot(df, aes(x = phy, y = micro)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Phylogenetic distance",
    y = "Microbiome distance"
  )

### ---- SPLICED FROM SCRIPT 2: species_colors -----------------------------
### (needed for cophylo tip/link coloring below)
species_colors <- c(
  # Cricetidae (blue shades)
  "Peromyscus melanotis"      = "#1f78b4",
  "Peromyscus boylii"         = "blue",
  "Peromyscus eremicus"       = "darkblue",
  "Microtus longicaudus"      = "#6baed6",
  "Reithrodontomys megalotis" = "#9ecae1",
  "Sigmodon ochrognathus"     = "#c6dbef",
  "Onychomys torridus"        = "#9ecae1",
  "Neotoma albigula"  = "#807dba",
  "Neotoma mexicana"  = "#6a5acd",
  
  # Heteromyidae (red shades)
  "Chaetodipus penicillatus"  = "#e41a1c",
  "Chaetodipus baileyi"       = "#fb6a4a",
  "Chaetodipus intermedius"   = "#fcae91",
  "Dipodomys merriami"        = "#a50f15",
  
  # Sciuridae (green shades)
  "Neotamias dorsalis"        = "#33a02c",
  "Ammospermophilus harrisii" = "#006d2c",
  
  # Soricidae (purple)
  "Sorex monticolus"          = "grey"
)
### -------------------------------------------------------------------------

### Cophylo --------------------------------------------------------------
assoc1 <- cbind(V1 = tree_exp$tip.label, V2 = tree_exp$tip.label)
species_colors1 <- species_colors[match(assoc1[, 1], names(species_colors))]

cophylo_obj <- cophylo(tree_exp, tree_avg, assoc1)

set.seed(9)
plot(
  cophylo_obj,
  type = "phylogram",
  link.type = "curved",
  link.lty = "solid",
  link.lwd = 2,
  link.col = species_colors1,
  ftype = "i",
  fsize = 0.7,
  gap = 100,
  mar = c(0.1, 0.1, 0.1, 0.1)
)


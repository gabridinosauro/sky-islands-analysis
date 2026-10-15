library(ape)
library(vegan)
library(dplyr)


#### BRMS model sky island project
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/")
bac.dist.sor = readRDS("bac.dist_sor.RDS")
bac.dist = readRDS("bac.dist.RDS")
sor_anaerobe = readRDS("sor_anaerobe.RDS")
sor_not = readRDS("sor_AEROTOL.RDS")
metadata = readRDS("/Users/gabri/Documents/final_codes_sky_islands/data/metadata_clean_1.RDS")

melted_BC <- reshape2::melt(as.matrix(bac.dist), varnames = c("Sample1", "Sample2"), value.name = "BC_Dissimilarity")
melted_sor = reshape2::melt(as.matrix(bac.dist.sor), varnames = c("Sample1", "Sample2"), value.name = "sor_Dissimilarity")
melted_an = reshape2::melt(as.matrix(sor_anaerobe), varnames = c("Sample1", "Sample2"), value.name = "sor_Dissimilarity_an")
melted_aert = reshape2::melt(as.matrix(sor_not), varnames = c("Sample1", "Sample2"), value.name = "sor_Dissimilarity_at")


unique(melted_BC$Sample1 == melted_sor$Sample1)
unique(melted_BC$Sample2 == melted_sor$Sample2) # Ok
unique(melted_BC$Sample1 == melted_an$Sample1)
unique(melted_BC$Sample2 == melted_an$Sample2) # Ok
unique(melted_BC$Sample1 == melted_aert$Sample1)
unique(melted_BC$Sample2 == melted_aert$Sample2) # Ok

melted_si = melted_sor  
melted_si$Sor_Similarity = 1 - melted_si$sor_Dissimilarity
melted_si$BC_Similarity =1 -  melted_BC$BC_Dissimilarity 
melted_si$sor_Similarity_an = 1- melted_an$sor_Dissimilarity_an 
melted_si$sor_Similarity_at = 1- melted_aert$sor_Dissimilarity_at 














names(melted_si)
#### remove half of the dataset
melted_si <- melted_si %>%
  dplyr::filter(Sample1 != Sample2) %>%
  mutate(pair_min = pmin(as.character(Sample1), as.character(Sample2)),
         pair_max = pmax(as.character(Sample1), as.character(Sample2))) %>%
  dplyr::distinct(pair_min, pair_max, .keep_all = TRUE) %>%
  dplyr::select(all_of(c("Sample1", "Sample2", "Sor_Similarity", "BC_Similarity","sor_Similarity_an","sor_Similarity_at")))

names(melted_si)


###### Add phylogentic distance
expanded_tree = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/expanded_tree.RDS")
expanded_tree$tip.label
host_phy_mat <- ape::cophenetic.phylo(expanded_tree)   





### Add phylogenetic distance
melted_si$phyl_dist <- host_phy_mat[
  cbind(as.character(melted_si$Sample1),
        as.character(melted_si$Sample2))
]

metadata$x_coord[1:10]
metadata$y_coord[1:10]
metadata$sample_id[1:10]

library(geosphere)

melted_si$geo_dist <- distm(
  cbind(metadata$x_coord, metadata$y_coord),
  fun = distHaversine
)[
  cbind(
    match(as.character(melted_si$Sample1), metadata$sample_id),
    match(as.character(melted_si$Sample2), metadata$sample_id)
  )
] / 1000

melted_si$elev_dist <- abs(
  metadata$elevation[match(as.character(melted_si$Sample1), metadata$sample_id)] -
    metadata$elevation[match(as.character(melted_si$Sample2), metadata$sample_id)]
)


mean(melted_si$phyl_dist)


# NDVI --------------
melted_si$NDVI_dist <- abs(
  metadata$mean_ndvi_2022_2025[match(as.character(melted_si$Sample1), metadata$sample_id)] -
    metadata$mean_ndvi_2022_2025[match(as.character(melted_si$Sample2), metadata$sample_id)]
)

# residualize NDVI distance on elevation distance
melted_si$NDVI_dist_resid <- residuals(lm(NDVI_dist ~ elev_dist, data = melted_si)) 

# Diet ---------------

nate_tree <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/dada2/Nates_tree.RDS")
nate_tree$edge.length <- nate_tree$edge.length[1:nrow(nate_tree$edge)]

phyl_mat <- ape::cophenetic.phylo(nate_tree)

diet_elt <- read.csv("/Users/gabri/Documents/final_codes_sky_islands/data/diet_elton.csv")

# 1. Phylogenetic distance matrix from tree, rename to match diet_elt format
rownames(phyl_mat) <- gsub("_", " ", rownames(phyl_mat))
colnames(phyl_mat) <- gsub("_", " ", colnames(phyl_mat))
rownames(phyl_mat)[rownames(phyl_mat) == "Tamias dorsalis"] <- "Neotamias dorsalis"
colnames(phyl_mat)[colnames(phyl_mat) == "Tamias dorsalis"] <- "Neotamias dorsalis"

# 2. Diet distance matrix
diet_mat <- as.matrix(diet_elt[, grep("^Diet\\.", colnames(diet_elt))])
rownames(diet_mat) <- diet_elt$tiplabel
diet_mat <- as.matrix(vegdist(diet_mat, method = "bray"))

# 3. Match species (drops Sorex monticolus, only in the tree)
common_sp <- intersect(rownames(phyl_mat), rownames(diet_mat))
length(common_sp)  # 15

phyl_mat <- phyl_mat[common_sp, common_sp]
diet_mat <- diet_mat[common_sp, common_sp]

# 4. Lower-triangle vectors, same species-pair order for both
phyl_vec <- phyl_mat[lower.tri(phyl_mat)]
diet_vec <- diet_mat[lower.tri(diet_mat)]

# 5. Regress diet on phylogeny, take residuals
lm_sp <- lm(diet_vec ~ phyl_vec)
summary(lm_sp)
diet_resid_vec <- residuals(lm_sp)
cor.test(diet_vec,phyl_vec)
# 6. Species-pair lookup table with BOTH raw and residualized diet distance,
#    for joining onto melted_si later
sp_pairs <- t(combn(common_sp, 2))
diet_lookup <- data.frame(
  sp1             = sp_pairs[, 1],
  sp2             = sp_pairs[, 2],
  diet_dist       = diet_vec,
  diet_dist_resid = diet_resid_vec
)

# 7. Map each sample to species, then join diet distances onto melted_si
melted_si$Species1 <- as.character(metadata$species[match(as.character(melted_si$Sample1), metadata$sample_id)])
melted_si$Species2 <- as.character(metadata$species[match(as.character(melted_si$Sample2), metadata$sample_id)])

key_melted <- paste(pmin(melted_si$Species1, melted_si$Species2),
                    pmax(melted_si$Species1, melted_si$Species2))
key_lookup <- paste(pmin(diet_lookup$sp1, diet_lookup$sp2),
                    pmax(diet_lookup$sp1, diet_lookup$sp2))

match_idx <- match(key_melted, key_lookup)

melted_si$diet_dist       <- diet_lookup$diet_dist[match_idx]
melted_si$diet_dist_resid <- diet_lookup$diet_dist_resid[match_idx]



melted_si$Species1 = metadata$species[match(as.character(melted_si$Sample1), metadata$sample_id)]
melted_si$Species2 = metadata$species[match(as.character(melted_si$Sample2), metadata$sample_id)]

melted_si$same_species = ifelse(melted_si$Species1 == melted_si$Species2, 1, 0)
melted_si = melted_si[melted_si$same_species ==0,]

melted_si$diet_resid_reg2 <- residuals(
  lm(diet_dist ~ phyl_dist, data = melted_si, na.action = na.exclude)
)

# should now be ~0 (that's the point)
cor.test(melted_si$diet_resid_reg2, melted_si$phyl_dist)


# standardize
cols <- c("phyl_dist", "geo_dist", "elev_dist", "NDVI_dist", "NDVI_dist_resid", "diet_dist","diet_dist", 
"diet_dist_resid", "diet_resid_reg2")
# Before scaling
sds <- sapply(melted_si[, cols], sd, na.rm = TRUE)
2 * sds



              melted_si <- melted_si %>%
  mutate(across(all_of(cols), ~ -((. - mean(., na.rm = TRUE)) /
                                    (2 * sd(., na.rm = TRUE)))))

getwd()




set.seed(123)  # optional, for reproducibility

sapply(melted_si[, cols], sd, na.rm = TRUE) * 2

melted_si$Mountain1 = metadata$Mountain[match(as.character(melted_si$Sample1), metadata$sample_id)]
melted_si$Mountain2 = metadata$Mountain[match(as.character(melted_si$Sample2), metadata$sample_id)]

metadata_sub <- metadata[sample(nrow(metadata), 100), ]

melted_si_sub <- melted_si[
  melted_si$Sample1 %in% metadata_sub$sample_id &
    melted_si$Sample2 %in% metadata_sub$sample_id,
]
melted_si$same_Mountain = ifelse(melted_si$Mountain1 == melted_si$Mountain2, 1, 0)

sds <- sapply(melted_si[, cols], sd, na.rm = TRUE)
2 * sds


melted_si$diet_dist
melted_si$phyl_dist
melted_si$elev_dist
melted_si$NDVI_dist
melted_si$NDVI_dist_resid
melted_si$diet_dist_resid

vars <- melted_si[, c("geo_dist", "phyl_dist", "elev_dist", "NDVI_dist", "NDVI_dist_resid", "diet_dist_resid", "diet_dist")]
cor_mat <- cor(vars, use = "pairwise.complete.obs")

full_names <- c("Geographic proximity", "Phylogenetic relateness", "Elevation similarity",
                "NDVI similarity", "NDVI similarity (resid.)", "Diet similarity (resid.)", "Diet similarity")

colnames(cor_mat) <- full_names
rownames(cor_mat) <- full_names

library(corrplot)

corrplot(cor_mat, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45,
         diag = FALSE,
         title = "Pairwise similarity (correlation)",
         mar = c(0, 0, 2, 0))

saveRDS(melted_si, "sky_islands_diet.RDS")
unique(melted_si$same_species)
melted_si$diet_dist_resid



### alpha and beat diversity at ASV level
##### Alpha and beta diversity 
#### libraries
library(vegan)
library(ggpubr)
library(ggplot2)
library(viridis)
source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")
#### setwd
setwd("/Users/gabri/Documents/final_codes_sky_islands/data")

#### read files in metadata
metadata = readRDS("metadata_clean_1.RDS")
ASVtable = readRDS("ASVtable.RDS")
dim(ASVtable)
Taxatable= readRDS("Taxtab.RDS")

### Check if samples in both tables are in the same order
unique(rownames(ASVtable) == metadata$sample_id) # OK
ASVtable = ASVtable[metadata$species != "Sorex monticolus", ] ### remove those samples from the ASvtable
metadata = metadata[metadata$species != "Sorex monticolus",]
unique(rownames(ASVtable) == metadata$sample_id) # OK


diet_elt <- read.csv("/Users/gabri/Documents/final_codes_sky_islands/data/diet_elton.csv")
# 2. Diet distance matrix
diet_mat <- as.matrix(diet_elt[, grep("^Diet\\.", colnames(diet_elt))])
rownames(diet_mat) <- diet_elt$tiplabel
diet_mat <- as.matrix(vegdist(diet_mat, method = "bray"))

# PCoA on diet Bray-Curtis distance matrix
pcoa_diet <- cmdscale(diet_mat, k = 2, eig = TRUE)

# Axis 1 scores per species
diet_axis1 <- pcoa_diet$points[, 1]
names(diet_axis1) <- rownames(diet_mat)
diet_axis1
unique(metadata$species)
metadata$diet_axis1 <- diet_axis1[as.character(metadata$species)]
metadata$diet_axis1



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
ASVtable <- ASVtable[, colSums(ASVtable) > 0]      # remove all-zero ASVs
ndvi_mod <- lm(mean_ndvi_2022_2025 ~ elevation, data = metadata)
metadata$ndvi_resid <- residuals(ndvi_mod)

###taxtab 
taxtab_fin = Taxatable[rownames(Taxatable) %in% colnames(ASVtable), ]


## ASVtable-----

### Calculate rarefied ASVtable
set.seed(12)
ASVtable_rar = rrarefy(ASVtable, 10000) ### rarefy the ASV table
#saveRDS(ASVtable_rar, "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")

### Calculate alpha diversity
metadata$Richness.rar<-specnumber(ASVtable_rar) # Counts ASV number
metadata$Shannon.rar<-vegan::diversity(ASVtable_rar, index="shannon") ## calculate Shannon H
##saveRDS(metadata,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/metadata.RDS")

m1 <- lm(Richness.rar ~ species + Mountain + elevation + ndvi_resid + diet_axis1, data = metadata)
anova(m1)
alias(m1)

metadata$diet_axis1
a = ggplot(metadata, aes(x = species, y = Richness.rar, color = Mountain)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.7) +
  theme_classic() +
  labs(y = "Number of ASVs", x = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
a


### check correlation to elevation
cor.test(metadata$Richness.rar,metadata$elevation, method = "pearson")


## check if quadratic model is better
m1 <- lm(Richness.rar ~ Mountain + elevation + ndvi_resid, data = metadata)
summary(m1)
m2 <- lm(Richness.rar ~ Mountain + poly(elevation, 2, raw = TRUE)+ ndvi_resid, data = metadata)
summary(m2)
AIC(m1, m2) ## M2 better



metadata$Species_ord = factor(metadata$species, levels = rev(levels(metadata$species)))
ggplot(metadata, aes(x = elevation, y = Richness.rar)) +
  facet_wrap(~factor(Mountain,
                     levels = c("Pinaleno", "Santa Catalina"),
                     labels = c("Pinaleño", "Santa Catalina")),
             ncol = 2)+
  geom_point(size = 2, alpha = 0.8 ,aes(color = Species_ord)) +
  scale_color_viridis_d(option = "turbo", begin = 0.05, end = 0.95) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              se = FALSE) +
  theme_classic() + ylab("Number of ASVs") + xlab("Elevation (m)")


## check if quadratic model is better
metadata_sc = metadata[metadata$Mountain == "Santa Catalina",]
metadata_pin = metadata[metadata$Mountain == "Pinaleno",]

### richness
m1 <- lm(Richness.rar ~  elevation, data = metadata_sc)
summary(m1)
m2 <- lm(Richness.rar ~  poly(elevation, 2, raw = TRUE), data = metadata_sc)
summary(m2)
AIC(m1, m2) ## M2 better for Santa catalina

### richness
m1 <- lm(Richness.rar ~  elevation, data = metadata_pin)
summary(m1)
m2 <- lm(Richness.rar ~  poly(elevation, 2, raw = TRUE), data = metadata_pin)
summary(m2)
AIC(m1, m2) ## M1 better for pinaleno


### Shannon
m1 <- lm(Shannon.rar ~  elevation, data = metadata_sc)
summary(m1)
m2 <- lm(Shannon.rar ~  poly(elevation, 2, raw = TRUE), data = metadata_sc)
summary(m2)
AIC(m1, m2) ## M2 better for Santa catalina

### Shannon
m1 <- lm(Shannon.rar ~  elevation, data = metadata_pin)
summary(m1)
m2 <- lm(Shannon.rar ~  poly(elevation, 2, raw = TRUE), data = metadata_pin)
summary(m2)
AIC(m1, m2) ## M1 better for pinaleno



### Altogheter
m1 <- lm(Shannon.rar ~  elevation, data = metadata)
summary(m1)
m2 <- lm(Shannon.rar ~  poly(elevation, 2, raw = TRUE), data = metadata)
summary(m2)
AIC(m1, m2) ## M2 better for pinaleno






AIC(m1, m2) ## M1 better

a = ggplot(metadata, aes(x = elevation, y = Richness.rar)) +
  #facet_wrap(~Mountain) +
  geom_point(size = 2, alpha = 0.8, aes(color = Species_ord)) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              se = FALSE) +
  labs(
    x = "Elevation (m)",
    y = "ASV Richness",
    color = "Species"
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  ) +  scale_color_viridis_d(option = "turbo", begin = 0.05, end = 0.95) 

b =  ggplot(metadata, aes(x = elevation, y = Shannon.rar)) +
  #facet_wrap(~Mountain) +
  geom_point(size = 2, alpha = 0.8, aes(color = Species_ord)) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              se = FALSE) +
  labs(
    x = "Elevation (m)",
    y = "Shannon H'",
    color = "Species"
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  ) +  scale_color_viridis_d(option = "turbo", begin = 0.05, end = 0.95) 

b

ggarrange(a,b, common.legend = TRUE, legend = "right")

## Beta diversity ----------------
set.seed(12)
#bac.dist<-vegdist(ASVtable_rar, method="bray") ### this calculates BC dissimilaries
#bac.dist.sor<-vegdist(ASVtable_rar, method="bray", binary  = TRUE ) ### this calculates Sorensen dissimilaries
#saveRDS(bac.dist,"bac.dist.RDS")
#saveRDS(bac.dist.sor,"bac.dist_sor.RDS")
bac.dist = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/bac.dist.RDS")
bac.dist.sor = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/bac.dist_sor.RDS")

set.seed(12)
bac.nmds<-metaMDS(bac.dist, k=2, try=999)
bac.nmds$stress #0.22
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

metadata$Axis01<-bac.nmds$points[,1]
metadata$Axis02<-bac.nmds$points[,2]



# install.packages("viridis") if needed
bac_nmds <- ggplot(metadata, aes(x = Axis01, y = Axis02)) +
  geom_point(aes(color = Species_ord, shape = Mountain), size = 2, alpha = 0.8) +
  stat_ellipse(aes(fill = Species_ord), geom = "polygon", alpha = 0.45, level = 0.68, show.legend = FALSE) +
  scale_color_viridis_d(option = "turbo", begin = 0.05, end = 0.95) +
  scale_fill_viridis_d(option  = "turbo", begin = 0.05, end = 0.95) +
  scale_shape_discrete(labels = c("Pinaleno" = "Pinale\u00f1o",
                                  "Santa Catalina" = "Santa Catalina")) +
  theme_classic() +
  labs(x = "NMDS1", y = "NMDS2", shape = "Mountain", color = "Species") +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  ) + annotate("text",
             x = Inf, y = Inf,
             label = stress,
             hjust = 1.1, vjust = 1.5,
             size = 3) + theme(
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.4, "cm")
          ) +
  guides(color = guide_legend(label.theme = element_text(size = 7, face = "italic")))

bac_nmds

perm = adonis2(bac.dist.sor~ species + Mountain + elevation +ndvi_resid +diet_axis1  , data = metadata, by = "margin")
perm


library(dplyr)
library(vegan)

results <- lapply(split(metadata, metadata$species), function(df) {
  df <- droplevels(df)
  
  # skip empty groups
  if (nrow(df) == 0) return(NULL)
  
  # need both mountains represented
  if (length(unique(df$Mountain)) < 2) return(NULL)
  
  # need enough samples overall
  if (nrow(df) < 10) return(NULL)
  
  # optional: require at least 2 samples per mountain
  tab <- table(df$Mountain)
  if (any(tab < 2)) return(NULL)
  
  dist_sub <- as.dist(as.matrix(bac.dist)[df$sample_id, df$sample_id])
  
  mod <- adonis2(dist_sub ~ Mountain, data = df, permutations = 999)
  
  data.frame(
    species = as.character(df$species[1]),
    n = nrow(df),
    n_mountains = length(unique(df$Mountain)),
    R2 = mod$R2[1],
    F = mod$F[1],
    p = mod$`Pr(>F)`[1]
  )
})

results_df <- bind_rows(results)
results_df
results_df <- do.call(rbind, results)
results_df
### plot to elevation 
table(metadata$mean_elev , metadata$species)



# mean elevation for each species
elev_df <- metadata %>%
  group_by(species) %>%
  summarise(mean_elev = unique(mean_elev)[1])

# merge with R² results
plot_df <- results_df %>%
  left_join(elev_df, by = "species")

plot_df

ggplot(plot_df, aes(x = mean_elev, y = R2)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  ggrepel::geom_text_repel(aes(label = species), size = 3) +
  labs(
    x = "Mean elevation (m)",
    y = expression(R^2),
    title = expression(paste("Relationship between mean elevation and ", R^2))
  ) +
  theme_bw()
cor.test(plot_df$R2, plot_df$mean_elev)

#anova(mod)
#saveRDS(perm,"perm.RDS")
#perm = readRDS("perm.RDS")
perm

ef <- envfit(bac.nmds, metadata["elevation"], permutations = 999)
ef
plot(bac.nmds); plot(ef)  




## Taxonomy ----
mod_tax = generate.tax.summary.modified(ASVtable_rar, taxtab_fin)
phyla = mod_tax$tax3
total_phyla = rowSums(phyla)
phyla_frac = total_phyla / sum(total_phyla)
sort(phyla_frac, decreasing = TRUE)


family_s = mod_tax$tax6
total_family_s = rowSums(family_s)
family_s_frac = total_family_s / sum(total_family_s)
sort(family_s_frac, decreasing = TRUE)

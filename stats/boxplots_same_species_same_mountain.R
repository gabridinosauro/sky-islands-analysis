### alpha and beat diversity at ASV level
##### Alpha and beta diversity 
#### libraries
library(vegan)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(viridis)
library(coin)
library(tidyverse)

source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")
#### setwd
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/")

#### read files in metadata
metadata = readRDS("metadata_clean.RDS")
ASVtable = readRDS("ASVtable.RDS")
Taxatable= readRDS("TAXtable_traits.RDS")



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
metadata = metadata[metadata$read_number >10000, ] ### remove those samples from the metadata
ASVtable = ASVtable[rowSums(ASVtable) >10000, ] ### remove those samples from the ASvtable

### Calculate rarefied ASVtable
set.seed(12)
ASVtable_rar = rrarefy(ASVtable, 10000) ### rarefy the ASV table
#saveRDS(ASVtable_rar, "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")



### same species different speices and combantion with mountain
sor_tot <- vegdist(ASVtable_rar, method = "bray", binary = TRUE)
# Convert to similarity
sor_tot <- 1 - as.matrix(sor_tot)

library(tidyverse)
sor_mat <- as.matrix(sor_tot)
sor_df <- as.data.frame(as.table(sor_tot))
colnames(sor_df) <- c("sample1", "sample2", "similarity")
# keep only upper triangle
sor_df <- sor_df[upper.tri(sor_mat), ]
any(sor_df$sample1 == sor_df$sample2)

meta <- metadata %>%
  select(sample_id, species, Mountain)


meta <- metadata %>%
  select(sample_id, species, Mountain)


sor_df <- sor_df %>%
  left_join(meta, by = c("sample1" = "sample_id")) %>%
  dplyr::rename(species1 = species, mountain1 = Mountain) %>%
  left_join(meta, by = c("sample2" = "sample_id")) %>%
  dplyr::rename(species2 = species, mountain2 = Mountain) %>%
  mutate(
    species_same  = ifelse(species1 == species2, "same_species", "different_species"),
    mountain_same = ifelse(mountain1 == mountain2, "same_mountain", "different_mountain")
  ) %>%
  mutate(
    species_same  = factor(species_same,  levels = c("same_species",  "different_species"),
                           labels = c("Same species", "Different species")),
    mountain_same = factor(mountain_same, levels = c("same_mountain", "different_mountain"),
                           labels = c("Same mountain", "Different mountain"))
  )

wilcox.test(as.numeric(sor_df$similarity)~ sor_df$species_same)

perm_test <- oneway_test(
  similarity ~ species_same,
  data = sor_df,
  distribution = approximate(B = 9999)
)

perm_test

summary_stats <- sor_df %>%
  group_by(species_same) %>%
  summarise(
    median = median(similarity, na.rm = TRUE),
    IQR = IQR(similarity, na.rm = TRUE),
    n = n()
  )

summary_stats

sor_df_sa_s = sor_df[sor_df$species_same == "Same species",]
wilcox.test(as.numeric(sor_df_sa_s$similarity)~ sor_df_sa_s$mountain_same)

perm_test <- oneway_test(
  similarity ~ mountain_same,
  data = sor_df_sa_s,
  distribution = approximate(B = 9999)
)

perm_test

sor_df_sa_s = sor_df[sor_df$species_same != "Same species",]
sor_df_sa_m = sor_df[sor_df$species_same == "Same species",]

wilcox.test(as.numeric(sor_df_sa_s$similarity)~ sor_df_sa_s$mountain_same)
library(coin)

perm_test <- oneway_test(
  similarity ~ mountain_same,
  data = sor_df_sa_s,
  distribution = approximate(B = 9999)
)

perm_test







a= ggplot(sor_df_sa_m, aes(x = mountain_same, y = similarity)) +
  #facet_wrap(~ species_same, strip.position = "top") +
  geom_violin(trim = TRUE, alpha = 0.6, color = NA, fill = "gray60") +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "gray20", fill = "white") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  ylab("Fraction of Shared ASVs") + xlab(NULL)# + coord_cartesian(ylim = c(0, 0.025)) 




a



b = ggplot(sor_df_sa_s, aes(x = mountain_same, y = similarity)) +
  #facet_wrap(~ species_same, strip.position = "top") +
  geom_violin(trim = TRUE, alpha = 0.6, color = NA, fill = "gray60") +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "gray20", fill = "white") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  ylab("Fraction of Shared ASVs") + xlab(NULL) + coord_cartesian(ylim = c(0, 0.1)) 

b


ggarrange(a,b)




#### Anaerboic -------------------------------
Taxatable1 <- Taxatable[match(colnames(ASVtable_rar), rownames(Taxatable)), , drop = FALSE]

# sanity check: all ASVs should match
all(rownames(Taxatable1) == colnames(ASVtable_rar))

keep_anaerobe <- !is.na(Taxatable1$anaerob) & Taxatable1$anaerob > 95
keep_aerobe   <- !is.na(Taxatable1$anaerob) & Taxatable1$anaerob < 5

fam_anaerobe <- ASVtable_rar[, keep_anaerobe, drop = FALSE]
fam_aerobe   <- ASVtable_rar[, keep_aerobe, drop = FALSE]

sum(is.na(fam_anaerobe))
sum(is.na(fam_aerobe))

bc_anaerobe <- vegdist(fam_anaerobe, method = "bray", binary = TRUE)
bc_aerobe   <- vegdist(fam_aerobe, method = "bray", binary = TRUE)


# Convert to similarity
sor_anaerobe <- 1 - as.matrix(bc_anaerobe)
sor_not <- 1 - as.matrix(bc_aerobe)
#saveRDS(sor_anaerobe, "sor_anaerobe.RDS")
#saveRDS(sor_not, "sor_AEROTOL.RDS")
getwd()
#sor_no_consensus <- 1 - as.matrix(bc_no_consensus)

### now I need between species
rownames(sor_anaerobe)
head(metadata)

sor_mat <- as.matrix(sor_anaerobe)
sor_df <- as.data.frame(as.table(sor_mat))
colnames(sor_df) <- c("sample1", "sample2", "similarity")
# keep only upper triangle
sor_df <- sor_df[upper.tri(sor_mat), ]
any(sor_df$sample1 == sor_df$sample2)

meta <- metadata %>%
  select(sample_id, species)

sor_df <- sor_df %>%
  left_join(meta, by = c("sample1" = "sample_id")) %>%
  dplyr::rename(species1 = species) %>%
  left_join(meta, by = c("sample2" = "sample_id")) %>%
  dplyr::rename(species2 = species) %>%
  mutate(
    species_same = ifelse(species1 == species2, "same_species", "different_species")
  )
sor_df$aerotellearnce = "obligate anaerobic"










# example with anaerobe (repeat for others if needed)
sor_mat <- as.matrix(sor_not)

sor_df1 <- as.data.frame(as.table(sor_mat))
colnames(sor_df1) <- c("sample1", "sample2", "similarity")

# keep only upper triangle
sor_df1 <- sor_df1[upper.tri(sor_mat), ]

any(sor_df1$sample1 == sor_df1$sample2)



meta <- metadata %>%
  select(sample_id, species)

sor_df1 <- sor_df1 %>%
  left_join(meta, by = c("sample1" = "sample_id")) %>%
  dplyr::rename(species1 = species) %>%
  left_join(meta, by = c("sample2" = "sample_id")) %>%
  dplyr::rename(species2 = species) %>%
  mutate(
    species_same = ifelse(species1 == species2, "same_species", "different_species")
  )

sor_df1$aerotellearnce = "aerotolerant"

sor_df_species = rbind(sor_df,sor_df1)
species = ggplot(sor_df_species, aes(x = aerotellearnce, y = similarity)) +
  geom_boxplot() + facet_wrap(~species_same) +
  theme_classic() +
  labs(
    x = "Aerotollerance",
    y = "Sørensen similarity",
    title = "Microbiome similarity within vs between species"
  )

species

perm_test <- function(df, group_var = "species_same", group_val, n_perm = 999) {
  d <- df[df[[group_var]] == group_val, ]
  obs <- mean(d$similarity[d$aerotellearnce == "aerotolerant"]) -
    mean(d$similarity[d$aerotellearnce == "obligate anaerobic"])
  perm <- replicate(n_perm, {
    sh <- sample(d$aerotellearnce)
    mean(d$similarity[sh == "aerotolerant"]) -
      mean(d$similarity[sh == "obligate anaerobic"])
  })
  list(observed = obs, p_value = mean(abs(perm) >= abs(obs)))
}
library(coin)
sor_df_species_same = sor_df_species[sor_df_species$species_same == "same_species",]
perm_test <- oneway_test(
  similarity ~ as.factor(aerotellearnce),
  data = sor_df_species_same,
  distribution = approximate(B = 9999)
)
perm_test


sor_df_species_same = sor_df_species[sor_df_species$species_same == "different_species",]
perm_test <- oneway_test(
  similarity ~ as.factor(aerotellearnce),
  data = sor_df_species_same,
  distribution = approximate(B = 9999)
)
perm_test

sor_df_species$species_same


##### MOuntain
sor_mat <- as.matrix(sor_anaerobe)

sor_df <- as.data.frame(as.table(sor_mat))
colnames(sor_df) <- c("sample1", "sample2", "similarity")

# keep only upper triangle
sor_df <- sor_df[upper.tri(sor_mat), ]

any(sor_df$sample1 == sor_df$sample2)

metadata$Mountain

meta <- metadata %>%
  select(sample_id, Mountain)

sor_df <- sor_df %>%
  left_join(meta, by = c("sample1" = "sample_id")) %>%
  dplyr::rename(Mountain1 = Mountain) %>%
  left_join(meta, by = c("sample2" = "sample_id")) %>%
  dplyr::rename(Mountain2 = Mountain) %>%
  mutate(
    Mountain_same = ifelse(Mountain1 == Mountain2, "same_Mountain", "different_Mountain")
  )

sor_df$aerotellearnce = "obligate anaerobic"


# example with anaerobe (repeat for others if needed)

sor_mat <- as.matrix(sor_not)

sor_df1 <- as.data.frame(as.table(sor_mat))
colnames(sor_df1) <- c("sample1", "sample2", "similarity")

# keep only upper triangle
sor_df1 <- sor_df1[upper.tri(sor_mat), ]

any(sor_df1$sample1 == sor_df1$sample2)



meta <- metadata %>%
  select(sample_id, Mountain)

sor_df1 <- sor_df1 %>%
  left_join(meta, by = c("sample1" = "sample_id")) %>%
  dplyr::rename(Mountain1 = Mountain) %>%
  left_join(meta, by = c("sample2" = "sample_id")) %>%
  dplyr::rename(Mountain2 = Mountain) %>%
  mutate(
    Mountain_same = ifelse(Mountain1 == Mountain2, "same_Mountain", "different_Mountain")
  )

sor_df1$aerotellearnce = "aerotolerant"

sor_df = rbind(sor_df,sor_df1)
Mountain = ggplot(sor_df, aes(x = aerotellearnce, y = similarity)) +
  geom_boxplot() + facet_wrap(~Mountain_same) +
  theme_classic() +
  labs(
    x = "Aerotollerance",
    y = "Sørensen similarity",
    title = "Microbiome similarity within vs between Mountain"
  )
Mountain

sor_df_mountain_same = sor_df[sor_df$Mountain_same == "same_Mountain",]
perm_test <- oneway_test(
  similarity ~ as.factor(aerotellearnce),
  data = sor_df_mountain_same,
  distribution = approximate(B = 9999)
)
perm_test


sor_df_mountain_same = sor_df[sor_df$Mountain_same == "different_Mountain",]
perm_test <- oneway_test(
  similarity ~ as.factor(aerotellearnce),
  data = sor_df_species_same,
  distribution = approximate(B = 9999)
)
perm_test

sor_df_species$species_same



ggarrange(species, Mountain, ncol = 1)

library(ggplot2)
library(ggpubr)

plot_aero <- function(df, facet_var) {
  if ("species_same" %in% names(df))
    df$species_same  = factor(df$species_same,  levels = c("same_species", "different_species"),
                              labels = c("Same species", "Different species"))
  if ("Mountain_same" %in% names(df))
    df$Mountain_same = factor(df$Mountain_same, levels = c("same_Mountain", "different_Mountain"),
                              labels = c("Same mountain", "Different mountain"))
  
  ggplot(df, aes(x = aerotellearnce, y = similarity, fill = aerotellearnce)) +
    geom_boxplot(width = 0.30, outlier.shape = NA, color = "gray20") +
    facet_wrap(as.formula(paste("~", facet_var)), ncol = 1, strip.position = "top") +
    scale_fill_manual(values = c("aerotolerant" = "#E07B54", "obligate anaerobic" = "#00BA38")) +
    coord_cartesian(ylim = c(0, quantile(df$similarity, 0.99))) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    ) +
    labs(x = NULL, y = "ASV sharing rate")
}
p1 <- plot_aero(sor_df_species, "species_same")
p1
p2 <- plot_aero(sor_df, "Mountain_same")
p2
ggarrange(p1, p2,  heights = c(1, 1), nrow = 1)






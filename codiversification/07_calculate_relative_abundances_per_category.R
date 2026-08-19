### things needed
#### Codiversification code HPC
### Load the libraries
library(vegan)
library(ape)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(tibble)
library(forcats)


`%nin%` = Negate(`%in%`)

### Load the files
### first load the ASV tree
### This will have to change once I loop through all my trees
#dir.create("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")
setwd("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")
tree_list = list.files("/Users/gabri/Desktop/mouse/SILVA_tree/subtrees", full.names = TRUE)
tree_host = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/expanded_tree.RDS")
ASVtab_rar = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")
Taxtab = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/Taxtab.RDS")
metadata_rar = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/metadata.RDS")
write.csv(metadata_rar,"/Users/gabri/Desktop/metadata.csv")


#### Calculate different fractions of co-diversification
results_final = readRDS("result_final.RDS") ### reads results final in
results_final$node_id = gsub("subtree_node" , ""  , results_final$node)
node_list = readRDS("nodes_list_1.RDS")
node_list$all_asvs<- lapply(node_list$descendants, names)
### Make lists of ASVs at different levels
results_final$p_cor = p.adjust(results_final$mantel_p_spear, method = "fdr")
results_final_sig = results_final[results_final$p_cor < 0.05,]
results_final_sig_weak  = results_final_sig[results_final_sig$mantel_stat_spear > 0.25,]
### all codiv ASVs (results final p_cor > 0.05)
node_list_weak_codiv = node_list[node_list$node_id %in% results_final_sig_weak$node_id,]
weak_codiv_asvs = unique(unlist(node_list_weak_codiv$all_asvs))

### medium codiv ASVs
results_final_sig_medium  = results_final_sig[results_final_sig$mantel_stat_spear > 0.50,]
### all codiv ASVs (results final p_cor > 0.05)
node_list_weak_medium = node_list[node_list$node_id %in% results_final_sig_medium$node_id,]
medium_codiv_asvs = unique(unlist(node_list_weak_medium$all_asvs))

### strong_codiv_ASVs
results_final_sig_strong  = results_final_sig[results_final_sig$mantel_stat_spear > 0.75,]
### all codiv ASVs (results final p_cor > 0.05)
node_list_weak_strong = node_list[node_list$node_id %in% results_final_sig_strong$node_id,]
strong_codiv_asvs = unique(unlist(node_list_weak_strong$all_asvs))

### 
weak_codiv_asvs = weak_codiv_asvs[weak_codiv_asvs %nin% medium_codiv_asvs]
medium_codiv_asvs = medium_codiv_asvs[medium_codiv_asvs %nin% strong_codiv_asvs]
length(weak_codiv_asvs) + length(medium_codiv_asvs) + length(strong_codiv_asvs)


#### ASVtab
ncol(ASVtab_rar)
colnames(ASVtab_rar)[1:10]
rownames(ASVtab_rar)
head(metadata_rar)


# 1. convert ASV table to long format
asv_long <- ASVtab_rar %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "abundance")


# 2. add codiversification info
asv_long <- asv_long %>%
  mutate(
    codiv_strength = case_when(
      ASV %in% strong_codiv_asvs ~ "strong",
      ASV %in% medium_codiv_asvs ~ "medium",
      ASV %in% weak_codiv_asvs ~ "weak",
      TRUE ~ "none"
    )
  )

# 3. sum abundance per sample × category
sample_category <- asv_long %>%
  group_by(sample_id, codiv_strength) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

# 4. compute relative abundance within each sample
sample_category <- sample_category %>%
  group_by(sample_id) %>%
  mutate(
    rel_abundance = abundance / sum(abundance, na.rm = TRUE)
  ) %>%
  ungroup()

# 3. add metadata (species)
sample_category <- sample_category %>%
  left_join(metadata_rar %>% select(sample_id, species),
            by = "sample_id")



ggplot(sample_category, aes(x = species, y = rel_abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~ codiv_strength) +
  labs(x = "Species", y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_category <- function(cat) {
  df <- sample_category %>%
    filter(codiv_strength == cat) %>%
    group_by(species) %>%
    summarise(
      med_rel_abundance = median(rel_abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(med_rel_abundance)) %>%
    mutate(species = factor(species, levels = species)) %>%
    left_join(
      sample_category %>% filter(codiv_strength == cat),
      by = "species"
    )
  
  ggplot(df, aes(x = species, y = rel_abundance)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    labs(x = "Species", y = "Relative abundance", title = cat) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_none   <- plot_category("none")
p_weak   <- plot_category("weak")
p_medium <- plot_category("medium")
p_strong <- plot_category("strong")

p_none
p_weak
p_medium
p_strong




strong_df <- sample_category %>%
  filter(codiv_strength == "strong")


ggplot(strong_df, aes(x = fct_reorder(species, rel_abundance, .fun = median, .desc = TRUE),
                      y = rel_abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





### Mean presence per specie per ASVs
asv_long <- ASVtab_rar %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "abundance") %>%
  left_join(metadata_rar %>% select(sample_id, species,Mountain),
            by = "sample_id")


asv_long <- asv_long %>%
  mutate(present = abundance > 0)
total_occupancy <- asv_long %>%
  group_by(ASV) %>%
  summarise(
    total_occupancy = sum(present, na.rm = TRUE),
    .groups = "drop"
  )%>%filter(total_occupancy > 0)%>%
  mutate(
    codiv_strength = case_when(
      ASV %in% strong_codiv_asvs ~ "strong",
      ASV %in% medium_codiv_asvs ~ "medium",
      ASV %in% weak_codiv_asvs ~ "weak",
      TRUE ~ "none"
    ))


results_final$mantel_stat_spear




species_occupancy <- asv_long %>%
  filter(present) %>%
  group_by(ASV) %>%
  summarise(
    species_occupancy = n_distinct(species),
    .groups = "drop"
  )%>%
  mutate(
    codiv_strength = case_when(
      ASV %in% strong_codiv_asvs ~ "strong",
      ASV %in% medium_codiv_asvs ~ "medium",
      ASV %in% weak_codiv_asvs ~ "weak",
      TRUE ~ "none"
    ))




species_occupancy_strong_none = species_occupancy[species_occupancy$codiv_strength %in% c("strong", "none"), ]
wilcox.test(species_occupancy ~ codiv_strength,
            data = species_occupancy_strong_none)

a = ggplot(species_occupancy_strong_none,
       aes(x = codiv_strength, y = species_occupancy)) +
    geom_jitter(width = 0.15, alpha = 0.01, size = 1) +
  geom_boxplot(outlier.shape = NA) +

  labs(x = "Codiversification Signal",
       y = "Species occupancy (per ASV)") +
  theme_classic()+ coord_cartesian(ylim = c(0, 8)) +scale_x_discrete(labels = c(
    "strong" = "Strong (rho > 0.75)",
    "none" = "None(rho < 0.25)"
  ))

a


total_occupancy_strong_none = total_occupancy[total_occupancy$codiv_strength %in% c("strong", "none"), ]

wilcox.test(total_occupancy ~ codiv_strength,
            data = total_occupancy_strong_none)

b = ggplot(total_occupancy_strong_none,
           aes(x = codiv_strength, y = total_occupancy)) +
  geom_jitter(width = 0.15, alpha = 0.01, size = 1) +
  geom_boxplot(outlier.shape = NA) +
  
  labs(x = "Codiversification Signal",
       y = "Sample occupancy (per ASV)") +
  theme_classic()+ coord_cartesian(ylim = c(0, 50))+scale_x_discrete(labels = c(
    "strong" = "Strong (rho > 0.75)",
    "none" = "None (rho < 0.25)"
  ))


b




total_occupancy_strong_none$divided = total_occupancy_strong_none$total_occupancy / species_occupancy_strong_none$species_occupancy
wilcox.test(divided ~ codiv_strength,
            data = total_occupancy_strong_none)
c = ggplot(total_occupancy_strong_none,
           aes(x = codiv_strength, y = divided)) +
  geom_jitter(width = 0.15, alpha = 0.01, size = 1) +
  geom_boxplot(outlier.shape = NA) +
  
  labs(x = "Codiversification Signal",
       y = "Total occupancy / Species occupancy (per ASV)") +
  theme_classic()+ coord_cartesian(ylim = c(0, 50))+scale_x_discrete(labels = c(
    "strong" = "Strong (rho > 0.75)",
    "none" = "None (rho < 0.25)"
  ))


c




ggpubr::ggarrange(a,b,c, nrow = 1)




### try 

Mountain_occupancy <- asv_long %>%
  filter(present) %>%
  group_by(ASV) %>%
  summarise(
    Mountain_occupancy = n_distinct(Mountain),
    .groups = "drop"
  )%>%
  mutate(
    codiv_strength = case_when(
      ASV %in% strong_codiv_asvs ~ "strong",
      ASV %in% medium_codiv_asvs ~ "medium",
      ASV %in% weak_codiv_asvs ~ "weak",
      TRUE ~ "none"
    ))




MOuntain_occupancy_strong_none = Mountain_occupancy[Mountain_occupancy$codiv_strength %in% c("strong", "none"), ]

wilcox.test(Mountain_occupancy ~ codiv_strength,
            data = MOuntain_occupancy_strong_none)

d = ggplot(MOuntain_occupancy_strong_none,
           aes(x = codiv_strength, y = Mountain_occupancy)) +
  geom_jitter(width = 0.15, alpha = 0.01, size = 1) +
  geom_boxplot(outlier.shape = NA) +
  
  labs(x = "Codiversification category",
       y = "Mountain occupancy (per ASV)") +
  theme_classic()+ coord_cartesian(ylim = c(0, 3))

d


### try correlations
#### read files in metadata
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/")
metadata = readRDS("metadata_clean.RDS")
ASVtable = readRDS("ASVtable.RDS")



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
ASVtable <- ASVtable[,colSums(ASVtable) != 0]
Taxatable1 = Taxatable[rownames(Taxatable) %in% colnames(ASVtable),]

rownames(ASVtable)
head(metadata)




# align metadata to ASVtable rows
meta2 <- metadata[match(rownames(ASVtable), metadata$sample_id), ]

# check alignment
stopifnot(all(rownames(ASVtable) == meta2$sample_id))

# count in how many host species each ASV is present
asv_species_count <- sapply(ASVtable, function(x) {
  length(unique(meta2$species[x > 0]))
})
### results nodelist 
head(results_final)
head(node_list)


node_list2 <- merge(
  node_list,
  results_final[, c("node_id", "mantel_stat_spear", "mantel_p_spear")],
  by = "node_id",
  all.x = TRUE
)

asv_codiv <- node_list2 %>%
  unnest(all_asvs) %>%
  group_by(all_asvs) %>%
  summarise(
    max_codiv = max(mantel_stat_spear, na.rm = TRUE),
    .groups = "drop"
  )
asv_codiv <- asv_codiv %>%
  filter(is.finite(max_codiv))

Taxatable_asvs <- Taxatable_asvs %>%
  tibble::rownames_to_column("ASV") %>%
  left_join(asv_codiv, by = c("ASV" = "all_asvs"))

Taxatable_asvs$genome
Taxatable_asvs$max_codiv = ifelse(Taxatable_asvs$max_codiv < 0, 0, Taxatable_asvs$max_codiv)
ggplot(Taxatable1, aes(x = species_count, y = max_codiv)) +
  geom_point(alpha = 0.2) +
  geom_smooth(
    method = "lm",
    #formula = y ~ x + I(x^2),
    se = FALSE
  ) +
  theme_classic() +
  labs(
    x = "Host Species Occupancy",
    y = "Codiversification signal per ASV"
  )


range01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
Taxatable1 = Taxatable_asvs
Taxatable1$aerob = 100 - Taxatable1$anaerob
Taxatable1$aerob_01 <- range01(Taxatable1$aerob)
Taxatable1$sporulation_01 <- range01(Taxatable1$sporulation)
Taxatable1$gram_01 <- range01(Taxatable1$gram)
Taxatable1$codiv <- range01(Taxatable1$max_codiv)
Taxatable1$species_count <- asv_species_count[Taxatable1$ASV]

Taxatable1$genome_z      <- as.numeric(scale(Taxatable1$genome))
Taxatable1$temp_z        <- as.numeric(scale(Taxatable1$Temperature))
Taxatable1$GC_z          <- as.numeric(scale(Taxatable1$GC))



model_std <- glm(
  cbind(species_count, 15 - species_count) ~
    aerob_01 +codiv+
    sporulation_01 +
    genome_z + I(genome_z^2) +
    temp_z + I(temp_z^2) +
    GC_z + I(GC_z^2) +
    gram_01,
  family = binomial,
  data = Taxatable1
)
summary(model_std)
plot(model_std)





#### model for 
summary_fam = as.data.frame(aggregate(species_count ~ Genus,
                        data = Taxatable1,
                        FUN = function(x)
                          c(mean = mean(x, na.rm = TRUE),
                            median = median(x, na.rm = TRUE),
                            sd = sd(x, na.rm = TRUE),
                            min = min(x, na.rm = TRUE),
                            max = max(x, na.rm = TRUE))))

summary_fam <- do.call(data.frame, summary_fam)

## remove all the junk from my species name 
con1 <- file("/Users/gabri/Documents/Sky_islands_project/gtdb_genus_summary.jsonl", open = "r")
data <- stream_in(con1)
close(con)
head(data)
str(data$summaries[[1]])


df <- lapply(seq_len(nrow(data)), function(i) {
  fam <- data$tax_name[i]
  s   <- data$summaries[[i]]
  
  # flatten the nested percentages dataframe
  s$pct_true  <- s$percentages$true
  s$pct_false <- s$percentages$false
  s$percentages <- NULL
  
  s$family <- fam
  s
}) |> bind_rows()


# aerobic -----------
df <- df |> select(family, everything())
df_ob_aer = df[df$name == "obligate anaerobic",]
df_ob_aer <- df_ob_aer |>
  mutate(verdict = case_when(
    grepl("^true",          majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false",         majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))



### 
summary_fam$anaerob = 100 - df_ob_aer$pct_false[match(summary_fam$Genus,df_ob_aer$family )]

summary_fam$anaerob = ifelse(is.na(summary_fam$anaerob), df_ob_aer$pct_true[match(summary_fam$Genus,df_ob_aer$family )],
                             summary_fam$anaerob)    

cor.test(summary_fam$anaerob,summary_fam$species_count.median)
plot(summary_fam$anaerob,summary_fam$species_count.median)



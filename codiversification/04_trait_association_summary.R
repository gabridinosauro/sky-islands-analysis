#### 04 - Family/phylum summary, GTDB trait associations, species comparison #
# Picks up from checkpoints saved in scripts 02 & 03.
# Part A: fraction of ASVs codiversifying per bacterial family
# Part B: correlate GTDB-derived life-history traits with codiversification
#         (obligate anaerobiosis, sporulation, genome size, growth temp, GC%, gram stain)
# Part C: phylum-level summary bar plot
# Part D: relative abundance of codiversifying taxa by host species
###############################################################################

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(colorspace)
library(jsonlite)
library(ggpubr)  # for ggarrange

setwd("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")

## ---- Load checkpoints -------------------------------------------------------
nodes_list        <- readRDS("nodes_list_1.RDS")
results_final     <- readRDS("results_final_padj.RDS")
tip_values_df      <- readRDS("tip_values_df.RDS")
Taxtab            <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/Taxtab.RDS")
ASVtab_rar        <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/ASVtable_rar.RDS")
metadata_rar      <- readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/analysis/data/metadata.RDS")

Taxtab1 <- Taxtab
Taxtab1$label <- rownames(Taxtab1)

## =============================================================================
## PART A - fraction codiversifying ASVs per family (top 50/100 by count)
## =============================================================================
merged_data <- merge(tip_values_df, Taxtab1, by = "label", all.x = TRUE)
merged_data$codiversifying <- !is.na(merged_data$mantel_final)
saveRDS(merged_data, "merged_data.RDS")

top50 <- merged_data %>%
  mutate(
    Family = str_squish(str_replace_all(as.character(Family), "[<>]", "")),
    Phylum = str_squish(str_replace_all(as.character(Phylum), "[<>]", "")),
    codiversifying = as.logical(codiversifying)
  ) %>%
  filter(!is.na(Family), Family != "", !str_detect(tolower(Family), "unknown")) %>%
  group_by(Family, Phylum) %>%
  summarise(
    n_ASVs = n(),
    n_codiv = sum(codiversifying, na.rm = TRUE),
    frac_codiv = n_codiv / n_ASVs,
    .groups = "drop"
  ) %>%
  arrange(desc(n_ASVs)) %>%
  slice_head(n = 100)

top50 <- top50[top50$n_ASVs > 10, ]

phyla <- unique(top50$Phylum)
phylum_colors <- setNames(
  darken(colorRampPalette(brewer.pal(12, "Set3"))(length(phyla)), 0.35),
  phyla
)
saveRDS(top50, "top50.RDS")
saveRDS(phylum_colors, "phylum_colors.RDS")

ggplot(top50, aes(x = reorder(Family, frac_codiv), y = frac_codiv, fill = Phylum)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = n_ASVs, color = Phylum), hjust = -0.05, size = 3) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = phylum_colors, name = "Phylum") +
  scale_color_manual(values = phylum_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                      expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Fraction codiversifying ASVs", title = "Co-diversification by bacteria family") +
  theme_classic() +
  theme(legend.position = "left", plot.margin = margin(5.5, 40, 5.5, 5.5))

## =============================================================================
## PART B - GTDB trait association plots (anaerobiosis, sporulation, genome
## size, growth temperature, GC%, gram stain) vs fraction codiversifying
## =============================================================================
con  <- file("/Users/gabri/Downloads/gtdb_family_summary.jsonl", open = "r")
data <- stream_in(con)
close(con)

df <- lapply(seq_len(nrow(data)), function(i) {
  fam <- data$tax_name[i]
  s   <- data$summaries[[i]]
  s$pct_true  <- s$percentages$true
  s$pct_false <- s$percentages$false
  s$percentages <- NULL
  s$family <- fam
  s
}) |> bind_rows() |> select(family, everything())
saveRDS(df, "gtdb_trait_long.RDS")

## -- obligate anaerobic
df_ob_aer <- df[df$name == "obligate anaerobic", ]
df_ob_aer <- df_ob_aer[df_ob_aer$family %in% top50$Family, ] %>%
  mutate(verdict = case_when(
    grepl("^true", majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false", majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))
df_ob_aer$pct_true <- ifelse(df_ob_aer$verdict == "false" & is.na(df_ob_aer$pct_true), 0, df_ob_aer$pct_true)

top50_1 <- top50 %>%
  left_join(df_ob_aer %>% select(family, majority_label, verdict, pct_true), by = c("Family" = "family"))
cor.test(top50_1$pct_true, top50$frac_codiv, method = "spearman")

a <- ggplot(top50_1, aes(x = pct_true, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains obligate anaerobic", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- sporulation
df_ob_aer <- df[df$name == "sporulation", ]
df_ob_aer <- df_ob_aer[df_ob_aer$family %in% top50$Family, ]
top50_1 <- top50 %>%
  left_join(df_ob_aer %>% select(family, majority_label, pct_true, pct_false), by = c("Family" = "family"))
top50_1$pct_true <- ifelse(is.na(top50_1$pct_true), 100 - top50_1$pct_false, top50_1$pct_true)
cor.test(top50_1$pct_true, top50$frac_codiv, method = "spearman")

b <- ggplot(top50_1, aes(x = pct_true, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains sporulating", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- genome size
df_gensize <- df[df$name %in% c("genome size", "estimated genome size"), ] %>%
  mutate(mean_val = as.numeric(gsub("Median: | bp", "", majority_label))) %>%
  group_by(family) %>%
  summarise(mean_val = mean(mean_val, na.rm = TRUE), .groups = "drop")
df_gensize <- df_gensize[df_gensize$family %in% top50$Family, ]
top50_1 <- top50 %>% left_join(df_gensize %>% select(family, mean_val), by = c("Family" = "family"))

c <- ggplot(top50_1, aes(x = mean_val, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  labs(x = "Mean genome size", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- growth temperature
df_gensize <- df[df$name %in% c("temperature growth"), ] %>%
  mutate(mean_val = as.numeric(gsub("Median: | Celsius", "", majority_label)))
df_gensize <- df_gensize[df_gensize$family %in% top50$Family, ]
top50_1 <- top50 %>% left_join(df_gensize %>% select(family, mean_val), by = c("Family" = "family"))

d <- ggplot(top50_1, aes(x = mean_val, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  labs(x = "Median optimal growth temperature (°C)", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- GC %
df_gensize <- df[df$name %in% c("GC percentage"), ] %>%
  mutate(mean_val = as.numeric(gsub("Median: | %", "", majority_label)))
df_gensize <- df_gensize[df_gensize$family %in% top50$Family, ]
top50_1 <- top50 %>% left_join(df_gensize %>% select(family, mean_val), by = c("Family" = "family"))

e <- ggplot(top50_1, aes(x = mean_val, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  labs(x = "GC %", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- gram positive
df_ob_aer <- df[df$name == "gram positive", ] %>%
  mutate(verdict = case_when(
    grepl("^true", majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false", majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))
top50_1 <- top50 %>%
  left_join(df_ob_aer %>% select(family, majority_label, pct_true), by = c("Family" = "family"))
top50_1$pct_true <- ifelse(is.na(top50_1$pct_true), 100 - top50_1$pct_false, top50_1$pct_true)
cor.test(top50_1$pct_true, top50$frac_codiv, method = "spearman")

f <- ggplot(top50_1, aes(x = pct_true, y = frac_codiv, color = Phylum)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains gram +", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

ggarrange(a, b, c, d, e, f, common.legend = TRUE, ncol = 2, nrow = 3)

## =============================================================================
## PART C - phylum-level codiversifying vs non-codiversifying summary
## =============================================================================
all_dominant <- unique(unlist(nodes_list$dominant_asvs, use.names = FALSE))
Taxtab_dom <- Taxtab1[rownames(Taxtab1) %in% all_dominant, ]

results_final_filt_new <- results_final[results_final$padj < 0.05, ]
results_final_filt_new <- results_final_filt_new[results_final_filt_new$mantel_stat_spear > 0.75, ]

codiv1 <- results_final_filt_new$node_name
codiv1 <- nodes_list[nodes_list$node_id %in% codiv1, ]
codiv1 <- unlist(codiv1$dominant_asvs)

Taxtab_dom$codiv <- ifelse(rownames(Taxtab_dom) %in% codiv1, "codiversifying", "non-codiversifying")
table(Taxtab_dom$codiv)

df_wide <- as.data.frame.matrix(table(Taxtab_dom$Phylum, Taxtab_dom$codiv))
df_wide$total <- rowSums(df_wide)
df_wide$fraction <- df_wide$codiversifying / df_wide$total
df_wide$Phylum <- rownames(df_wide)
df_wide$Phylum <- factor(df_wide$Phylum, levels = df_wide$Phylum[order(df_wide$fraction, decreasing = TRUE)])

ggplot(df_wide, aes(x = Phylum, y = fraction)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = codiversifying), vjust = -1.5, size = 3) +
  geom_text(aes(label = total), vjust = -0.3, color = "orange", size = 3) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = "Phylum", y = "Fraction ASVs codiversifying") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## =============================================================================
## PART D - relative abundance of codiversifying taxa, by host species
## =============================================================================
asv_long <- ASVtab_rar %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "abundance")

asv_long <- asv_long %>%
  left_join(Taxtab_dom %>% rownames_to_column("ASV") %>% select(ASV, codiv), by = "ASV") %>%
  left_join(metadata_rar %>% select(sample_id, species), by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(rel_abundance = abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup()

codiv_only <- asv_long %>% filter(codiv == "codiversifying")

codiv_sum <- codiv_only %>%
  group_by(sample_id, species) %>%
  summarise(
    total_codiv_abundance = sum(abundance, na.rm = TRUE),
    total_codiv_rel_abundance = sum(rel_abundance, na.rm = TRUE),
    .groups = "drop"
  )

species_order <- codiv_sum %>%
  group_by(species) %>%
  summarise(mean_rel = mean(total_codiv_rel_abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_rel))

codiv_sum$species <- factor(codiv_sum$species, levels = species_order$species)

ggplot(codiv_sum, aes(x = species, y = total_codiv_rel_abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(x = "Host Species", y = "Relative abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

cat("Done with trait association & summary plots.\n")

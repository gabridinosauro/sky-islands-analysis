#### 06 - GTDB trait permutation tests ########################################
# Picks up top50.RDS, phylum_colors.RDS and gtdb_trait_long.RDS from script 04.
# For each GTDB life-history trait (obligate anaerobiosis, sporulation, genome
# size, growth temperature, GC%, gram stain), tests whether codiversifying
# families are enriched/depleted for that trait using a permutation approach:
# the codiversification weight (n_codiv, or frac_codiv for sporulation - see
# note below) is randomly reshuffled across families 10,000 times to build a
# null distribution, which the observed weighted statistic is compared against.
#
# TWO THINGS CARRIED OVER FROM YOUR ORIGINAL EXPLORATORY CODE - WORTH A LOOK:
#   1. The sporulation test weights by `frac_codiv`, while every other trait
#      test weights by `n_codiv` (raw codiversifying ASV count). Kept as-is
#      below since I don't know if that was intentional; flagging it here.
#   2. Your original script's phylum-level codiversifying summary bar chart
#      (further downstream in your draft) filtered nodes on
#      `mantel_stat_pear` (Pearson r), whereas script 04's Part C used
#      `mantel_stat_spear` (Spearman r) for the same figure. That figure is
#      already produced in script 04 so it isn't repeated here - just decide
#      which correlation type should define "codiversifying" for it.
###############################################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggpubr)

setwd("/Users/gabri/Desktop/mouse/SILVA_tree/treeees")

## ---- Load checkpoints (script 04) -------------------------------------------
top50         <- readRDS("top50.RDS")
phylum_colors <- readRDS("phylum_colors.RDS")
df            <- readRDS("gtdb_trait_long.RDS")

## =============================================================================
## Rebuild the six trait-vs-codiversification scatterplots
## (these are the "main" panel shown above each null-distribution plot below)
## =============================================================================

## -- obligate anaerobic
df_anaer <- df[df$name == "obligate anaerobic", ]
df_anaer <- df_anaer[df_anaer$family %in% top50$Family, ] %>%
  mutate(verdict = case_when(
    grepl("^true", majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false", majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))
df_anaer$pct_true <- ifelse(df_anaer$verdict == "false" & is.na(df_anaer$pct_true), 0, df_anaer$pct_true)

top50_anaer <- top50 %>%
  left_join(df_anaer %>% select(family, majority_label, verdict, pct_true), by = c("Family" = "family"))

p_anaerobic <- ggplot(top50_anaer, aes(x = pct_true, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains obligate anaerobic", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- sporulation
df_spor <- df[df$name == "sporulation", ]
df_spor <- df_spor[df_spor$family %in% top50$Family, ]
top50_spor <- top50 %>%
  left_join(df_spor %>% select(family, majority_label, pct_true, pct_false), by = c("Family" = "family"))
top50_spor$pct_true <- ifelse(is.na(top50_spor$pct_true), 100 - top50_spor$pct_false, top50_spor$pct_true)

p_sporulation <- ggplot(top50_spor, aes(x = pct_true, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains sporulating", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- genome size
df_gsize <- df[df$name %in% c("genome size", "estimated genome size"), ] %>%
  mutate(mean_val = as.numeric(gsub("Median: | bp", "", majority_label))) %>%
  group_by(family) %>%
  summarise(mean_val = mean(mean_val, na.rm = TRUE), .groups = "drop")
df_gsize <- df_gsize[df_gsize$family %in% top50$Family, ]
top50_gsize <- top50 %>% left_join(df_gsize %>% select(family, mean_val), by = c("Family" = "family"))

p_genome_size <- ggplot(top50_gsize, aes(x = mean_val, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  labs(x = "Mean genome size", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- growth temperature
df_temp <- df[df$name %in% c("temperature growth"), ] %>%
  mutate(mean_val = as.numeric(gsub("Median: | Celsius", "", majority_label)))
df_temp <- df_temp[df_temp$family %in% top50$Family, ]
top50_temp <- top50 %>% left_join(df_temp %>% select(family, mean_val), by = c("Family" = "family"))

p_temperature <- ggplot(top50_temp, aes(x = mean_val, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  labs(x = "Median optimal growth temperature (°C)", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- GC%
df_gc <- df[df$name %in% c("GC percentage"), ] %>%
  mutate(mean_val = as.numeric(gsub("Median: | %", "", majority_label)))
df_gc <- df_gc[df_gc$family %in% top50$Family, ]
top50_gc <- top50 %>% left_join(df_gc %>% select(family, mean_val), by = c("Family" = "family"))

p_gc <- ggplot(top50_gc, aes(x = mean_val, y = frac_codiv)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  labs(x = "GC %", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

## -- gram stain
df_gram <- df[df$name == "gram positive", ] %>%
  mutate(verdict = case_when(
    grepl("^true", majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false", majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))
top50_gram <- top50 %>%
  left_join(df_gram %>% select(family, majority_label, pct_true, pct_false), by = c("Family" = "family"))
top50_gram$pct_true <- ifelse(is.na(top50_gram$pct_true), 100 - top50_gram$pct_false, top50_gram$pct_true)

p_gram <- ggplot(top50_gram, aes(x = pct_true, y = frac_codiv, color = Phylum)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains gram +", y = "Fraction co-diversifying ASVs", color = "Phylum") +
  theme_classic() + scale_color_manual(values = phylum_colors)

ggarrange(p_anaerobic, p_sporulation, p_genome_size, p_temperature, p_gc, p_gram,
          common.legend = TRUE, ncol = 2, nrow = 3)

## =============================================================================
## Permutation tests
## =============================================================================
set.seed(42)
B <- 10000

## -- obligate anaerobic: weighted mean %anaerobic (upper-tail test)
dd_anaer  <- top50_anaer[!is.na(top50_anaer$pct_true) & !is.na(top50_anaer$frac_codiv), ]
obs_anaer  <- weighted.mean(dd_anaer$pct_true, w = dd_anaer$n_codiv)
null_anaer <- replicate(B, weighted.mean(dd_anaer$pct_true, w = sample(dd_anaer$n_codiv)))
p_anaer_upper <- mean(null_anaer >= obs_anaer)
cat("Obligate anaerobic: observed =", round(obs_anaer, 2),
    " null mean =", round(mean(null_anaer), 2),
    " one-sided p (upper tail) =", p_anaer_upper, "\n")

p_null_anaer <- ggplot(data.frame(null_anaer), aes(null_anaer)) +
  geom_density(fill = "grey80", colour = NA) +
  geom_vline(xintercept = obs_anaer, linetype = 2, colour = "#c0392b", linewidth = 0.8) +
  geom_vline(xintercept = mean(null_anaer), linetype = 3, colour = "grey40") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains obligate anaerobic", y = "Null") +
  theme_classic()

(p_anaerobic + labs(x = NULL) + theme(axis.text.x = element_blank())) / p_null_anaer +
  plot_layout(heights = c(4, 1))

## -- sporulation: weighted mean %sporulating
## NOTE: weighted by frac_codiv, not n_codiv (see note at top of script)
d_spor <- top50_spor[!is.na(top50_spor$pct_true) & !is.na(top50_spor$frac_codiv), ]
obs_spor  <- weighted.mean(d_spor$pct_true, w = d_spor$frac_codiv)
null_spor <- replicate(B, weighted.mean(sample(d_spor$pct_true), w = d_spor$frac_codiv))
p_spor_lower <- (sum(null_spor <= obs_spor) + 1) / (B + 1)
p_spor_two   <- (sum(abs(null_spor - mean(null_spor)) >= abs(obs_spor - mean(null_spor))) + 1) / (B + 1)
cat("Sporulation: observed =", round(obs_spor, 2),
    " null mean =", round(mean(null_spor), 2),
    " one-sided p (lower tail) =", p_spor_lower,
    " two-sided p =", p_spor_two, "\n")

p_null_spor <- ggplot(data.frame(null_spor), aes(null_spor)) +
  geom_density(fill = "grey80", colour = NA) +
  geom_vline(xintercept = obs_spor, linetype = 2, colour = "#c0392b", linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "% strains sporulating", y = "Null") +
  theme_classic()

(p_sporulation + labs(x = NULL) + theme(axis.text.x = element_blank())) / p_null_spor +
  plot_layout(heights = c(4, 1))

## -- genome size: weighted SD among codiversifying families
wvar <- function(x, w) { w <- w / sum(w); mu <- sum(w * x); sum(w * (x - mu)^2) }

d_gsize <- top50_gsize[!is.na(top50_gsize$mean_val), ]
obs_gsize_sd  <- sqrt(wvar(d_gsize$mean_val, d_gsize$n_codiv))
null_gsize_sd <- replicate(B, sqrt(wvar(d_gsize$mean_val, sample(d_gsize$n_codiv))))
p_gsize <- (sum(null_gsize_sd <= obs_gsize_sd) + 1) / (B + 1)
cat("Genome size: observed weighted SD =", round(obs_gsize_sd, 0),
    " null mean SD =", round(mean(null_gsize_sd), 0),
    " p =", p_gsize, "\n")

## -- growth temperature: weighted mean distance from 37°C (host body temp)
dd_temp <- top50_temp[!is.na(top50_temp$mean_val) & !is.na(top50_temp$frac_codiv), ]
target  <- 37
obs_temp  <- weighted.mean(abs(dd_temp$mean_val - target), w = dd_temp$n_codiv)
null_temp <- replicate(B, weighted.mean(abs(dd_temp$mean_val - target), w = sample(dd_temp$n_codiv)))
p_temp <- (sum(null_temp <= obs_temp) + 1) / (B + 1)
cat("Growth temperature: observed mean |T-37| =", round(obs_temp, 2),
    " null mean =", round(mean(null_temp), 2),
    " p =", p_temp, "\n")

## -- GC%: weighted mean GC among codiversifying families (two-sided)
dd_gc <- top50_gc[!is.na(top50_gc$mean_val) & !is.na(top50_gc$frac_codiv), ]
obs_gc  <- weighted.mean(dd_gc$mean_val, w = dd_gc$n_codiv)
null_gc <- replicate(B, weighted.mean(dd_gc$mean_val, w = sample(dd_gc$n_codiv)))
p_gc_two <- (sum(abs(null_gc - mean(null_gc)) >= abs(obs_gc - mean(null_gc))) + 1) / (B + 1)
cat("GC%: observed =", round(obs_gc, 2),
    " null mean =", round(mean(null_gc), 2),
    " two-sided p =", p_gc_two, "\n")

## -- gram stain: weighted mean %gram-positive (two-sided)
dd_gram <- top50_gram[!is.na(top50_gram$pct_true) & !is.na(top50_gram$frac_codiv), ]
obs_gram  <- weighted.mean(dd_gram$pct_true, w = dd_gram$n_codiv)
null_gram <- replicate(B, weighted.mean(dd_gram$pct_true, w = sample(dd_gram$n_codiv)))
p_gram_two <- (sum(abs(null_gram - mean(null_gram)) >= abs(obs_gram - mean(null_gram))) + 1) / (B + 1)
cat("Gram stain: observed =", round(obs_gram, 2),
    " null mean =", round(mean(null_gram), 2),
    " two-sided p =", p_gram_two, "\n")

## ---- Summary table + checkpoint ---------------------------------------------
perm_summary <- data.frame(
  trait = c("obligate_anaerobic", "sporulation", "genome_size_sd",
            "growth_temperature", "GC_percent", "gram_positive"),
  observed  = c(obs_anaer, obs_spor, obs_gsize_sd, obs_temp, obs_gc, obs_gram),
  null_mean = c(mean(null_anaer), mean(null_spor), mean(null_gsize_sd),
                mean(null_temp), mean(null_gc), mean(null_gram)),
  p_value   = c(p_anaer_upper, p_spor_two, p_gsize, p_temp, p_gc_two, p_gram_two)
)
print(perm_summary)

saveRDS(perm_summary, "trait_permutation_results.RDS")
cat("Done. Saved trait_permutation_results.RDS\n")

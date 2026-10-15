### Code analysis BRMS sky island #
library(brms)
library(dplyr)
library(marginaleffects)
mm = readRDS("/Users/gabri/Documents/final_codes_sky_islands/data/Sor_Similarity_NDVI_diet_res2.RDS")
summary(mm)
performance::check_collinearity(mm)

post <- brms::as_draws_df(mm)

# look at posterior correlations
cor(post[, c("b_phyl_dist","b_geo_dist", "b_elev_dist", "b_NDVI_dist_resid", "b_diet_resid_reg2")])


model2 = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/dada2/Sor_Similarity_NDVI.RDS")

summary(model2)
#ordbetareg::pp_check_ordbeta(model2) 
post <- as_draws_df(model2)
# look at posterior correlations
cor(post[, c("b_phyl_dist", "b_geo_dist", "b_elev_dist","b_NDVI_dist_resid")]) #good
### check VIF
performance::check_collinearity(model2)

data_clean <- as.data.frame(model2$data)
#data_clean <- data_clean[, !sapply(data_clean, is.list)]
#sl <- marginaleffects::avg_slopes(model2, newdata = data_clean, ndraws = 500)
saveRDS(sl, "sl.RDS")

sl$term <- factor(sl$term, 
                   levels = c("phyl_dist", "elev_dist", "NDVI_dist_resid","geo_dist"),
                   labels = c("Phylogenetic relatedness","Elevation similarity","Vegetation cover similarity", "Geographic proximity"))

sl$category = "whole community"


ggplot(sl, aes(x = estimate, y = term)) +
  geom_point(
    size = 1,
    position = position_dodge(width = 0.4)
  ) +
  geom_errorbarh(
    aes(xmin = conf.low - 0.00015, xmax = conf.high + 0.00015), ## 0.000015 added for visualization otherwsie peopel think there is no error vars (sorry)
    height = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  labs(
    x = "Δ ASV-sharing rate (per 1 SD change in predictor)",
    y = NULL,
    title = "Cross-species ASV sharing"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 9, hjust = 1),plot.title = element_text(size = 10)
  )


###  aerolorelant
model_an = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/dada2/Sor_Similarity_at.RDS")

summary(model_an)
#ordbetareg::pp_check_ordbeta(model2) 
post <- as_draws_df(model_an)
# look at posterior correlations
cor(post[, c("b_phyl_dist", "b_geo_dist", "b_elev_dist","b_NDVI_dist_resid")]) #good
### check VIF
performance::check_collinearity(model_an)

data_clean <- as.data.frame(model_an$data)
#data_clean <- data_clean[, !sapply(data_clean, is.list)]
#sl1 <- avg_slopes(model_an, newdata = data_clean, ndraws = 500)
#saveRDS(sl1, "/Users/gabri/Documents/final_codes_sky_islands/data/sl1_at.RDS")
sl1 = readRDS("/Users/gabri/Documents/final_codes_sky_islands/data/sl1_at.RDS")
sl1$term <- factor(sl1$term, 
                  levels = c("phyl_dist", "elev_dist", "NDVI_dist_resid","geo_dist"),
                  labels = c("Phylogenetic relatedness","Elevation similarity","Vegetation cover similarity", "Geographic proximity"))


sl1[, 3:5] <- -sl1[, 3:5]
sl1$category = "Aerotolerant"

### anaerobic -----
model_an = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/dada2/Sor_Similarity_an.RDS")

summary(model_an)
#ordbetareg::pp_check_ordbeta(model2) 
post <- as_draws_df(model_an)
# look at posterior correlations
cor(post[, c("b_phyl_dist", "b_geo_dist", "b_elev_dist")]) #good
### check VIF
performance::check_collinearity(model_an)

data_clean <- as.data.frame(model_an$data)
#data_clean <- data_clean[, !sapply(data_clean, is.list)]
#sl2 <- avg_slopes(model_an, newdata = data_clean, ndraws = 500)
#saveRDS(sl2, "/Users/gabri/Documents/final_codes_sky_islands/data/sl2_an.RDS")
sl2 = readRDS("/Users/gabri/Documents/final_codes_sky_islands/data/sl2_an.RDS")

sl2$term <- factor(sl2$term, 
                   levels = c("phyl_dist", "elev_dist", "NDVI_dist_resid","geo_dist"),
                   labels = c("Phylogenetic relatedness","Elevation similarity","Vegetation cover similarity", "Geographic proximity"))

sl2[, 3:5] <- -sl2[, 3:5]
sl2$category = "Anaerobic"




sl_fin = rbind(sl1,sl2)

pd <- position_dodge(width = 0.5)

ggplot(sl_fin, aes(x = estimate, y = term, color = category)) +
  geom_errorbar(aes(xmin = conf.low - 0.0005, xmax = conf.high+ 0.0005),
                width = 0.5, orientation = "y", position = pd) +
  geom_point(size = 1.5, position = pd) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Δ ASV-sharing rate (per 1 SD change in predictor)", y = NULL,
       title = "Cross-species ASV sharing") +
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust = 1))


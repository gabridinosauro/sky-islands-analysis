# code to run ordbetareg on my computer
cat("R script started/n")
flush.console()
library(brms)
library(ordbetareg)
library(marginaleffects)
cat("brms loaded/n")
flush.console()

library(cmdstanr)
cat("cmdstanr loaded/n")
flush.console()


Sys.setenv(OMP_NUM_THREADS = 6)
cat("Sys.setenv set to 6 threads/n")
flush.console()

setwd("/xdisk/laurenpetrullo/schiro/sky_island/brms_models")
cat("setwd done/n")
flush.console()

# read_RDS
md_sky_islands = readRDS("sky_islands.RDS")
md_sky_islands = md_sky_islands[md_sky_islands$same_species ==0,]

model1<-ordbetareg(sor_Similarity_an~1+
                   phyl_dist+geo_dist+elev_dist+NDVI_dist_resid+
                   (1|mm(Sample1,Sample2))+(1|mm(Species1,Species2)),
                   data = md_sky_islands,
                   #family= "Beta",
                   warmup = 1000, iter = 3000,
                   threads = threading(4),
                   cores = 4,
                   chains = 4,
                   init=0)

saveRDS(model1,"Sor_Similarity_an.RDS")

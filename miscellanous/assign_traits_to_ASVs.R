### Assign traits to ASVs
### first genus
### aerotolerant ----------
library(jsonlite)
Taxatable= readRDS("Taxtab.RDS")



### first try species---------
### remove all the junk from my species name 
Taxatable$Species_clean =  gsub("_", " ", sub("\\(.*", "", Taxatable$Species))


con <- file("/Users/gabri/Downloads/gtdb_species_summary.jsonl", open = "r")
data <- stream_in(con)
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
Taxatable$anaerob = 100 - df_ob_aer$pct_false[match(Taxatable$Species_clean,df_ob_aer$family )]
Taxatable$anaerob = ifelse(is.na(Taxatable$anaerob), df_ob_aer$pct_true[match(Taxatable$Species_clean,df_ob_aer$family )],
                           Taxatable$anaerob)    




###sporulation ------
# move family to front
df <- df |> select(family, everything())
df_ob_aer = df[df$name == "sporulation",]
df_ob_aer <- df_ob_aer |>
  mutate(verdict = case_when(
    grepl("^true",          majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false",         majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))




Taxatable$sporulation = df_ob_aer$pct_true[match(Taxatable$Species_clean,df_ob_aer$family )]
Taxatable$sporulation = ifelse(is.na(Taxatable$sporulation), 100 - df_ob_aer$pct_false[match(Taxatable$Species_clean,df_ob_aer$family )],
                           Taxatable$anaerob)




###gram ------
# move family to front
df <- df |> select(family, everything())
df_ob_aer = df[df$name == "gram positive",]
df_ob_aer <- df_ob_aer |>
  mutate(verdict = case_when(
    grepl("^true",          majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false",         majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))





Taxatable$gram = df_ob_aer$pct_true[match(Taxatable$Species_clean,df_ob_aer$family )]
Taxatable$gram = ifelse(is.na(Taxatable$gram), 100 - df_ob_aer$pct_false[match(Taxatable$Species_clean,df_ob_aer$family )],
                               Taxatable$gram)





### genome size ---
df_spor = df[df$name == "genome size",]
df_spor$genome_median <- as.numeric(gsub(".*Median: ([0-9.]+) bp.*", "\\1", df_spor$majority_label))
### 
Taxatable$genome = df_spor$genome_median[match(Taxatable$Species_clean,df_spor$family )]



### GC content ----
df_gensize <- df[df$name %in% c("GC percentage"), ] |>
  mutate(mean_val = as.numeric(gsub("Median: | %", "", majority_label)))
Taxatable$GC = df_gensize$mean_val[match(Taxatable$Species_clean,df_gensize$family )]



### OPtimal growth temperature ----
df_gensize <- df[df$name %in% c("temperature growth"), ] |>
  mutate(mean_val = as.numeric(gsub("Median: | Celsius", "", majority_label)))
Taxatable$Temperature = df_gensize$mean_val[match(Taxatable$Species_clean,df_gensize$family )]










### genus ---------
## remove all the junk from my species name 
con1 <- file("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/dada2/gtdb_genus_summary.jsonl", open = "r")
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
Taxatable$anaerob = ifelse(is.na(Taxatable$anaerob),100 - df_ob_aer$pct_false[match(Taxatable$Genus,df_ob_aer$family )],
                           Taxatable$anaerob)


Taxatable$anaerob = ifelse(is.na(Taxatable$anaerob), df_ob_aer$pct_true[match(Taxatable$Genus,df_ob_aer$family )],
                           Taxatable$anaerob)    




###sporulation ------
# move family to front
df <- df |> select(family, everything())
df_ob_aer = df[df$name == "sporulation",]
df_ob_aer <- df_ob_aer |>
  mutate(verdict = case_when(
    grepl("^true",          majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false",         majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))




Taxatable$sporulation = ifelse(is.na(Taxatable$sporulation),df_ob_aer$pct_true[match(Taxatable$Genus,df_ob_aer$family )],
                               Taxatable$sporulation)
Taxatable$sporulation = ifelse(is.na(Taxatable$sporulation), 100 - df_ob_aer$pct_false[match(Taxatable$Genus,df_ob_aer$family )],
                               Taxatable$sporulation)




###gram ------
# move family to front
df <- df |> select(family, everything())
df_ob_aer = df[df$name == "gram positive",]
df_ob_aer <- df_ob_aer |>
  mutate(verdict = case_when(
    grepl("^true",          majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false",         majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))





Taxatable$gram = ifelse(is.na(Taxatable$gram),df_ob_aer$pct_true[match(Taxatable$Genus,df_ob_aer$family )],
                        Taxatable$gram)
Taxatable$gram = ifelse(is.na(Taxatable$gram), 100 - df_ob_aer$pct_false[match(Taxatable$Genus,df_ob_aer$family )],
                        Taxatable$gram)





### genome size
df_spor = df[df$name == "genome size",]
df_spor$genome_median <- as.numeric(gsub(".*Median: ([0-9.]+) bp.*", "\\1", df_spor$majority_label))
### 
Taxatable$genome = ifelse(is.na(Taxatable$genome),
                          df_spor$genome_median[match(Taxatable$Genus,df_spor$family )],
                          Taxatable$genome)



### GC content
df_gensize <- df[df$name %in% c("GC percentage"), ] |>
  mutate(mean_val = as.numeric(gsub("Median: | %", "", majority_label)))
Taxatable$GC = ifelse(is.na(Taxatable$GC),
  df_gensize$mean_val[match(Taxatable$Genus,df_gensize$family )],
  Taxatable$GC)



### OPtimal growth temperature
df_gensize <- df[df$name %in% c("temperature growth"), ] |>
  mutate(mean_val = as.numeric(gsub("Median: | Celsius", "", majority_label)))
Taxatable$Temperature = df_gensize$mean_val[match(Taxatable$Genus,df_gensize$family )]


Taxatable$Temperature = ifelse(is.na(Taxatable$Temperature),
                               df_gensize$mean_val[match(Taxatable$Genus,df_gensize$family )],
                               Taxatable$Temperature)














### family ---------
## remove all the junk from my species name 
con1 <- file("/Users/gabri/Documents/Sky_islands_project/gtdb_family_summary.jsonl", open = "r")
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
Taxatable$anaerob = ifelse(is.na(Taxatable$anaerob),100 - df_ob_aer$pct_false[match(Taxatable$Family,df_ob_aer$family )],
                           Taxatable$anaerob)


Taxatable$anaerob = ifelse(is.na(Taxatable$anaerob), df_ob_aer$pct_true[match(Taxatable$Family,df_ob_aer$family )],
                           Taxatable$anaerob)    




###sporulation ------
# move family to front
df <- df |> select(family, everything())
df_ob_aer = df[df$name == "sporulation",]
df_ob_aer <- df_ob_aer |>
  mutate(verdict = case_when(
    grepl("^true",          majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false",         majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))




Taxatable$sporulation = ifelse(is.na(Taxatable$sporulation),df_ob_aer$pct_true[match(Taxatable$Family,df_ob_aer$family )],
                               Taxatable$sporulation)
Taxatable$sporulation = ifelse(is.na(Taxatable$sporulation), 100 - df_ob_aer$pct_false[match(Taxatable$Family,df_ob_aer$family )],
                               Taxatable$sporulation)




###gram ------
# move family to front
df <- df |> select(family, everything())
df_ob_aer = df[df$name == "gram positive",]
df_ob_aer <- df_ob_aer |>
  mutate(verdict = case_when(
    grepl("^true",          majority_label, ignore.case = TRUE) ~ "true",
    grepl("^false",         majority_label, ignore.case = TRUE) ~ "false",
    grepl("No robust majority", majority_label, ignore.case = TRUE) ~ "no_robust_majority"
  ))





Taxatable$gram = ifelse(is.na(Taxatable$gram),df_ob_aer$pct_true[match(Taxatable$Family,df_ob_aer$family )],
                        Taxatable$gram)
Taxatable$gram = ifelse(is.na(Taxatable$gram), 100 - df_ob_aer$pct_false[match(Taxatable$Family,df_ob_aer$family )],
                        Taxatable$gram)





### genome size
df_spor = df[df$name == "genome size",]
df_spor$genome_median <- as.numeric(gsub(".*Median: ([0-9.]+) bp.*", "\\1", df_spor$majority_label))
### 
Taxatable$genome = ifelse(is.na(Taxatable$genome),
                          df_spor$genome_median[match(Taxatable$Family,df_spor$family )],
                          Taxatable$genome)



### GC content
df_gensize <- df[df$name %in% c("GC percentage"), ] |>
  mutate(mean_val = as.numeric(gsub("Median: | %", "", majority_label)))
Taxatable$GC = ifelse(is.na(Taxatable$GC),
                      df_gensize$mean_val[match(Taxatable$Family,df_gensize$family )],
                      Taxatable$GC)



### OPtimal growth temperature
df_gensize <- df[df$name %in% c("temperature growth"), ] |>
  mutate(mean_val = as.numeric(gsub("Median: | Celsius", "", majority_label)))
Taxatable$Temperature = df_gensize$mean_val[match(Taxatable$Family,df_gensize$family )]


Taxatable$Temperature = ifelse(is.na(Taxatable$Temperature),
                               df_gensize$mean_val[match(Taxatable$Genus,df_gensize$family )],
                               Taxatable$Temperature)
saveRDS(Taxatable,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/TAXtable_traits.RDS")

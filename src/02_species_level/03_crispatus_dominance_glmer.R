# This script evaluates the association between maternal and daughter vaginal 
# microbiome composition, specifically focusing on Lactobacillus crispatus dominance.

# Dependences: R 4.4.0, tidyverse, tidytacos

library(tidytacos)
library(tidyverse)
library(lme4)
library(performance)

# Load tt object derived from "01_isala_species_level_profiling"

isala <- 
  read_tidytacos(
    "results/isala/species_level_analyses/tt_isala_shotgun_processed/") %>%
  set_rank_names(
    c("kingdom","phylum","class","order","family","genus","species")) %>%
  add_taxon_name(include_species = TRUE) 

isala <-
  isala %>% 
  mutate_taxa(taxon_name=str_replace(taxon_name, "s__", ""),
              taxon_name=str_replace(taxon_name, "_", " "))  

# Clean metadata
samples_inf <- 
  samples(isala) %>%
  select(code, family_number, sample_type, dominant_sp, dominant_sp_2)

samples_inf <- 
  samples_inf %>%
  mutate(
    crispatus_dom = ifelse(
      dominant_sp == "Lactobacillus crispatus" |
        dominant_sp_2 == "Lactobacillus crispatus",
      1, 0
    )
  )

mother_status <- 
  samples_inf %>%
  filter(sample_type == "M") %>%
  select(family_number, mother_crispatus = crispatus_dom)

daughters_status <- 
  samples_inf %>%
  filter(sample_type %in% c("D", "D1", "D2", "D3", "D4")) %>%
  left_join(mother_status, by = "family_number")

# Fit generalized linear mixed model

model_glmer <- 
  glmer(
  crispatus_dom ~ mother_crispatus + (1 | family_number),
  data = daughters_status,
  family = binomial
)

# Inspect results

summary(model_glmer)

# p = 0.00597 **

# Odds ratios for fixed effects
exp(fixef(model_glmer))

# OR: 20.9999946

# Compute R2 (Nakagawa & Schielzeth)

r2_results <- r2(model_glmer)

r2_results

# Sensitivity analysis with glm

model_glm <- 
  glm(
  crispatus_dom ~ mother_crispatus,
  data = daughters_status,
  family = binomial
)

summary(model_glm)
exp(coef(model_glm))

# Similar p and odd ratio -> The association is robust to the inclusion of 
# family‑level random effects.
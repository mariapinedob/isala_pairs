# This script plots the taxonomic profiles of the samples from the non-EU cohort
# (France et al. 2022) as well as extracts data of important variables 
# (i.e., relative abundance of species of interest).

# Load required libraries
library(tidyverse)
library(tidytacos)
library(janitor)
library(ggplot2)

# Define path of output directory
dout <- "results/non_eu_cohort/species_level_analyses/"

# Create output directory if it doesn't exist
if (!dir.exists(dout)) {
  dir.create(dout, recursive = TRUE)
}

# Load tidytacos object
tt <- 
  read_tidytacos("data/tt_non_eu_cohort/") %>%
  mutate_samples(
    sample = str_replace_all(sample, "_1.clean_1.fastq.gz", ""))

# Load metadata and join
codes <- read_csv("data/non_eu_metadata.csv") %>%
  rename(sample = srr)

tt <- 
  tt %>%
  add_metadata(codes) %>%
  mutate_samples(
    code_type = str_extract(participant_code, "([MD]|D\\d{1,2}|GM|GD)$"))

# Clean taxon names
tt <- 
  tt %>%
  add_taxon_name(include_species = TRUE) %>%
  mutate_taxa(
    taxon_name = str_replace(taxon_name, "s__", ""),
    taxon_name = str_replace(taxon_name, "_", " ")
  )

# Retrieve relative abundance matrix
rel_abundance <- 
  rel_abundance_matrix(tt, sample_name = sample, taxon_name = taxon_name)

# Retrieve information of species present (>1%)
species_present <- 
  as.data.frame(rel_abundance) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "taxon_name", values_to = "abundance") %>%
  group_by(sample) %>%
  filter(abundance > 0.01)

ravel_samples <- samples(tt)
species_present <- species_present %>% left_join(ravel_samples, by = "sample")

# Retrieve dominant species information
dominance_thr <- 0.29
dominant_species <- as.data.frame(rel_abundance) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "taxon_name", values_to = "abundance") %>%
  group_by(sample) %>%
  summarise(
    max_sp1 = taxon_name[which.max(abundance)],
    max_sp1_rel_abundance = round(abundance[which.max(abundance)], 3),
    max_sp2 = taxon_name[order(abundance, decreasing = TRUE)[2]],
    max_sp2_rel_abundance = round(abundance[order(abundance, decreasing = TRUE)[2]], 3),
    dominant_sp = ifelse(max_sp1_rel_abundance >= dominance_thr, max_sp1, "No dominance"),
    dominant_sp_2 = ifelse(max_sp2_rel_abundance >= dominance_thr, max_sp2, "No co-dominance"),
    .groups = 'drop'
  )

# Relative abundance of species of interest

species_interest <- 
  rel_abundance %>%
  as.data.frame() %>%
  mutate(sample = rownames(rel_abundance)) %>%
  select("sample", "Lactobacillus crispatus", "Lactobacillus iners", 
         "Bifidobacterium leopoldii", "Bifidobacterium vaginale") %>%
  as_tibble()

# Integrate metadata into tt
tt <- 
  tt %>%
  add_metadata(dominant_species) %>%
  add_metadata(species_interest)

samples <- 
  samples(tt) %>%
  clean_names() %>%
  select(participant_code, bifidobacterium_leopoldii, bifidobacterium_vaginale,
         lactobacillus_crispatus, lactobacillus_iners)

write_tsv(samples, "data/ravel_sylph.tsv")

# Save updated tidytacos object
write_tidytacos(tt, dout = paste0(dout, "tidytacos_france_et_al"))

# Visualization settings
tt <- 
  tt %>%
  mutate_samples(
    code_type = factor(
      code_type, levels = c("M", "D", "D1", "D2", "D3", "D4")))

# Color palette
isala_colors <- 
  c("#BEA5D5", "#FF9896", "#D62728", "#8C564B", "#FF7F0E",  
                   "#FFBB78", "#1F77B4", "#2CA02C", "#AEC7E8", "#FAF3DD")

# Plot taxonomic profiles
plot_taxa <- 
  tt %>%
  tacoplot_stack(x = code_type, order_by = "code_type", n = 10) +
  facet_wrap(~family_code, nrow = 3, scales = "free_x") +
  theme_bw() +
  xlab("Participants") +
  ylab("Relative abundance (%)") +
  guides(fill = guide_legend(title = "Species")) +
  scale_fill_manual(values = isala_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "italic"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.spacing.x = unit(1.5, "pt"),
    panel.grid = element_blank()
  )

plot_taxa

# Save plot
ggsave(
  paste0(dout, "taxonomic_profiles.svg"),
  plot_taxa,
  units = "cm", width = 30, height = 14, dpi = 300
)

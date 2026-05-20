# This script plots the taxonomic profiles of the samples as well as extracts
# data of important variables (i.e., relative abundance of species of interest).

# Dependences: R 4.4.0, tidyverse, tidytacos

library(tidyverse)
library(tidytacos)
library(ggplot2)
library(janitor)

# Define path of output directory
dout <- "results/isala/species_level_analyses/"

# Create output directory if it doesn't exist
if (!dir.exists(dout)) {
  dir.create(dout, recursive = TRUE)
}

# Load tidytacos file
isala <- 
  read_tidytacos("data/tt_isala_shotgun_sylph/")

# Load metadata file with quality info from workflow 01 and merge
metadata <- 
  read_csv("data/metadata_quality_info.csv") %>%
  rename(sample = bioSampleId)

isala <- 
  isala %>% add_metadata(metadata)

# Filter high-quality samples
isala <- 
  isala %>% 
  filter_samples(high_quality_shotgun != FALSE)

# Save tt object with high-quality samples
write_tidytacos(isala, "data/tt_isala_shotgun_sylph_hq")

# Clean taxon names
isala <- 
  isala %>%
  set_rank_names(
    c("kingdom","phylum","class","order","family","genus","species")) %>%
  add_taxon_name(include_species = TRUE) %>% 
  mutate_taxa(taxon_name=str_replace(taxon_name, "s__", ""),
              taxon_name=str_replace(taxon_name, "_", " ")) 

# Retrieve relative abundance matrix
rel_abundance <- 
  rel_abundance_matrix(
    isala, sample_name = sample, taxon_name = taxon_name
    )

# All species

species_present <- 
  as.data.frame(rel_abundance) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "taxon_name", values_to = "abundance") %>%
  group_by(sample) %>%
  filter(abundance > 0) 

# Dominant species

dominance_thr <- 0.29

dominant_species <- 
  as.data.frame(rel_abundance) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "taxon_name", values_to = "abundance") %>%
  group_by(sample) %>%
  summarise(
    max_sp1 = taxon_name[which.max(abundance)],
    max_sp1_rel_abundance = round(abundance[which.max(abundance)], 3),
    max_sp2 = taxon_name[order(abundance, decreasing = TRUE)[2]],
    max_sp2_rel_abundance = round(abundance[order(abundance, decreasing = TRUE)[2]], 3),
    dominant_sp = ifelse(max_sp1_rel_abundance >= dominance_thr, 
                         max_sp1,
                         "No dominance"),
    dominant_sp_2 = ifelse(max_sp2_rel_abundance >= dominance_thr, 
                           max_sp2,
                           "No co-dominance"),
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

# Integrate in tt of Isala

isala <-
  isala %>%
  add_metadata(dominant_species) 

isala <- 
  isala %>%
  add_metadata(species_interest)  

# save metadata for supplementary tables

suppl_table1 <-
  samples(isala) %>%
  clean_names() %>%
  select(code, lactobacillus_crispatus, lactobacillus_iners,
         bifidobacterium_vaginale, bifidobacterium_leopoldii)

#write_csv(suppl_table1, "data/suppl_table1.csv")

# Extract other metadata for visualization
isala <- 
  isala %>%
  mutate_samples(
    family_number = str_extract(code, "F[0-9]{1,2}"),
    sample_type = str_extract(code, "(GM|M|GD|D\\d{0,2})"),
    code_type = factor(sample_type, 
                       levels = c("GM", "M", "GD", "D", "D1", "D2", "D3", "D4")
                       )
  )

# Save updated tidytacos object
write_tidytacos(isala, dout = paste0(dout, "tt_isala_shotgun_processed"))

order_plot <- c("GM", "M", "GD", "D", "D1", "D2", "D3", "D4")
  
isala_colors <- 
  c("#BEA5D5", "#8C564B", "#FF9896", "#D62728", "#FF7F0E",
                   "#FFBB78", "#1F77B4", "#AEC7E8", "#2CA02C", "#FAF3DD")

# Plot taxonomic profiles
isala %>%
  tacoplot_stack(x = code_type, order_by = "code_type", n = 12) +
  facet_wrap(~family_number, nrow = 3, scales = "free_x") +
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

# Save plot
ggsave(
  paste0(dout, "taxonomic_profiles_isala.svg"),
  units = "cm", width = 28, height = 14, dpi = 300
)

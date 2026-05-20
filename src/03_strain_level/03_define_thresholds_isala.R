# This script processes inStrain compare outputs derived from metagenomes in the
# Isala cohort, converts popANI values into population SNVs per megabase 
# (popSPM) and assesses their distributions across comparisons within and 
# between families.

# Dependences: R 4.4.0, tidyverse

library(tidyverse)
library(readr)
library(readxl)
library(janitor)

# Load directories
dir <- "results/isala/strain_level_analyses/inStrain_compare/"
dout <- file.path(dir, "results", "processed")

dir.create(dout, recursive = TRUE, showWarnings = FALSE)

# Load genome-wide comparison output (single or multiple files)
files <- 
  list.files(dir, pattern = "genomeWide_compare\\.tsv$", full.names = TRUE)

all_genome_info <- 
  map_dfr(
  files,
  ~ read_tsv(.x, show_col_types = FALSE)
)

# Extract sample structure from instrain naming and annotate comparisons
all_genome_info <- 
  all_genome_info %>%
  mutate(
    sample1 = str_remove(name1, "^ref_genomes_vs_"),
    sample2 = str_remove(name2, "^ref_genomes_vs_"),
    
    participant1 = str_extract(sample1, "^F\\d+"),
    participant2 = str_extract(sample2, "^F\\d+"),
    
    participant1_fam = str_extract(participant1, "F\\d{1,2}"),
    participant2_fam = str_extract(participant2, "F\\d{1,2}")
  ) %>%
  distinct(sample1, sample2, .keep_all = TRUE) %>%
  mutate(
    comparison_type = case_when(
      participant1_fam == participant2_fam ~ "same family",
      TRUE ~ "different family"
    ),
    spm = (1 - popANI) * 1e6
  )

# Prepare for visualization
plot_data <- 
  all_genome_info %>%
  mutate(
    genome = factor(genome, levels = c(
      "L_crispatus.fna",
      "L_iners.fna",
      "G_vaginalis.fna",
      "G_leopoldii.fna"
    )),
    genome_label = fct_recode(
      genome,
      "Lactobacillus crispatus" = "L_crispatus.fna",
      "Lactobacillus iners" = "L_iners.fna",
      "Gardnerella vaginalis" = "G_vaginalis.fna",
      "Gardnerella leopoldii" = "G_leopoldii.fna"
    ),
    pair_label = case_when(
      comparison_type == "same family" ~ "Same family",
      comparison_type == "different family" ~ "Different family"
    )
  )

# Visualization
p <- 
  plot_data %>%
  ggplot(aes(x = spm, fill = pair_label, color = pair_label)) +
  geom_histogram(alpha = 0.5, bins = 200) +
  facet_wrap(~ genome_label, scales = "free", nrow = 2) +
  scale_fill_manual(values = c("#8EB8A7", "#F5D491")) +
  scale_color_manual(values = c("#8EB8A7", "#F5D491")) +
  scale_x_log10() +
  geom_vline(xintercept = 100, linetype = "dashed", color = "blue") +
  labs(
    x = "popSPM",
    y = "Frequency",
    fill = "",
    color = ""
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    strip.text = element_text(face = "bold.italic", size = 12),
    legend.position = "bottom",
    panel.spacing.x = unit(5, "lines")
  )

p

ggsave(
  file.path(dout, "distribution_metagenomes_only.svg"),
  plot = p,
  width = 30,
  height = 15,
  units = "cm",
  dpi = 600
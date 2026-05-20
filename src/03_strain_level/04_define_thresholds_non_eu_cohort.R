# This script processes inStrain compare outputs from metagenomes of the
# non-European cohort (France et al. 2022), converts popANI values into 
# population SNVs per megabase (popSPM), and evaluates their distributions 
# across comparisons within and between families.

# Dependences: R 4.4.0, tidyverse

library(tidyverse)

# Load directories

dir <- "results/non_eu_cohort/strain_level_analyses/inStrain_compare/"
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

# Load metadata linking SRR to participant
codes <- 
  read_csv("data/non_eu_metadata.csv") %>%
  select(srr, participant_code)

# Extract sample IDs and map to participants
all_genome_info <- 
  all_genome_info %>%
  mutate(
    sample1 = str_extract(name1, "(?<=ref_genomes_vs_)[^_]+"),
    sample2 = str_extract(name2, "(?<=ref_genomes_vs_)[^_]+")
  ) %>%
  left_join(codes, by = c("sample1" = "srr")) %>%
  rename(participant1 = participant_code) %>%
  left_join(codes, by = c("sample2" = "srr")) %>%
  rename(participant2 = participant_code)

# Annotate comparisons
all_genome_info <- all_genome_info %>%
  distinct(sample1, sample2, .keep_all = TRUE) %>%
  mutate(
    participant1_fam = str_extract(participant1, "^[A-Z]\\d+"),
    participant2_fam = str_extract(participant2, "^[A-Z]\\d+"),
    comparison_type = case_when(
      participant1 == participant2 ~ "same individual",
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
  file.path(dout, "distribution_all_ravel.svg"),
  plot = p,
  width = 30,
  height = 15,
  units = "cm",
  dpi = 600
)
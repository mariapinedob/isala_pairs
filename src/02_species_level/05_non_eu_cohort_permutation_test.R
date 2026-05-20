# This script calculates Bray-Curtis dissimilarity between mother-daughter pairs
# and performs a permutation test on metagenomes from the cohort of France et al
# (2022)

# Load required libraries
library(tidyverse)
library(tidytacos)
library(readxl)

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
codes <- 
  read_csv("data/non_eu_metadata.csv") %>%
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

# Calculate Bray-Curtis dissimilarity index for real mother-daughter pairs
betas <- 
  betas(tt, unique = FALSE) %>%
  filter(
    sample_type_pt_1 == "mother", 
    family_code_1 == family_code_2, beta != 0) %>%
  select(
    sample_1, participant_code_1, sample_2, participant_code_2, beta)

observed_beta_values <- betas %>% select(beta)

###################################################
########## Permutation matrix #####################
###################################################

num_permutations <- 10000
permuted_matrix <- matrix(nrow = 25, ncol = num_permutations)

for (i in 1:num_permutations) {
  mothers <- tt %>% filter_samples(tt$samples$sample_type_pt == "mother")
  daughters <- tt %>% filter_samples(tt$samples$sample_type_pt == "daughter")
  
  # Randomize mother-daughter pairing
  daughters$samples$mother_pt2 <- sample(daughters$samples$mother_pt, replace = FALSE)
  names(daughters$samples$mother_pt2) <- daughters$samples$mother_pt
  
  # Merge and compute betas
  tt2 <- merge_tidytacos(mothers, daughters, taxon_identifier = taxon_name)
  permutedbeta <- betas(tt2, unique = FALSE) %>%
    filter(sample_1 == mother_pt2_2) %>%
    select(beta)
  
  permuted_matrix[, i] <- permutedbeta$beta
}

# Save permutation matrix
write.csv(permuted_matrix, "data/non_eu_permutation_matrix.csv", row.names = FALSE)

###################################################
################ Visualization #####################
###################################################

# Load permutation matrix
fin_matrix <- "data/non_eu_permutation_matrix.csv"
permuted_matrix <- read.csv(fin_matrix)

# Density plot of mean beta values
tibble(mean_beta = colMeans(permuted_matrix)) %>%
  ggplot(aes(x = mean_beta)) +
  geom_density() +
  geom_rug(data = tibble(mean_beta_real = mean(observed_beta_values[[1]])), aes(x = mean_beta_real)) +
  xlim(c(0, 1))

# Summary stats
mean_random <- mean(colMeans(permuted_matrix))
mean_real <- mean(observed_beta_values[[1]])
p_value <- sum(colMeans(permuted_matrix) < mean_real) / num_permutations

# Similarity values
obs_similarity <- 1 - observed_beta_values
similarity_matrix <- 1 - permuted_matrix

# Violin plot comparison
violin_plot_comparison <- 
  tibble(mean_beta = colMeans(similarity_matrix), group = "Random pairs") %>%
  ggplot(aes(x = group, y = mean_beta)) +
  geom_violin(fill = "#FCD0BE", color = "black", adjust = 3.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_segment(data = tibble(mean_beta_real = mean(obs_similarity[[1]]), group = "Real pairs"),
               aes(x = group, xend = group, y = mean_beta_real, yend = mean_beta_real),
               color = "#00719A", size = 1) +
  geom_point(data = tibble(mean_beta = obs_similarity[[1]], group = "Real pairs"),
             aes(x = group, y = mean_beta), color = "#00719A", size = 1.2, alpha = 0.5) +
  geom_segment(data = tibble(mean_beta_real = mean(obs_similarity[[1]])),
               aes(x = 1.5, xend = 2, y = mean_beta_real, yend = mean_beta_real),
               color = "#00719A", size = 1) +
  geom_text(data = tibble(mean_beta_real = mean(obs_similarity[[1]]), group = "Real pairs"),
            aes(x = group, y = mean_beta_real, label = round(mean_beta_real, 3)),
            color = "black", size = 4, hjust = -0.3) +
  theme_bw() +
  xlab("") +
  ylab("Mean Bray-Curtis Similarity Index (10K permutations)") +
  scale_x_discrete(limits = c("Real pairs", "Random pairs"),
                   labels = c("Random pairs" = "Average similarity \nof randomized pairs",
                              "Real pairs" = "Observed similarity \nof true pairs")) +
  scale_y_continuous(breaks = seq(0.0, 1.0, by = 0.2)) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.ticks.x = element_blank())

print(violin_plot_comparison)

# Save violin plot
ggsave(paste0(dout, "/permutation_test.svg"), violin_plot_comparison,
       units = "cm", width = 12, height = 14, dpi = 600)

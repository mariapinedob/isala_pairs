# This script generates a permutation matrix and performs a comparison
# between true and randomized mother-daughters pairs from Isala. 

# Dependences: R 4.4.0, tidyverse, tidytacos

library(tidyverse)
library(tidytacos)

# Define path of output directory
dout <- "results/isala/species_level_analyses/"

# Create output directory if it doesn't exist
if (!dir.exists(dout)) {
  dir.create(dout, recursive = TRUE)
}

# Load tacos object: Only samples with high quality
isala <- 
  read_tidytacos("data/tt_isala_shotgun_sylph_hq/") %>%
  add_taxon_name()

# Create variable "sample_type_pt" (Sample type in permutation test)
isala <-
  isala %>%
  mutate_samples(
    family_number = str_extract(code, "F\\d{1,2}"),
    sample_type   = str_extract(code, "GM|M|GD|D\\d{0,2}"),
    code_type     = factor(sample_type,
                           levels = c("GM", "M", "GD", "D", "D1", "D2", "D3", "D4")),
    sample_type_pt = if_else(code_type == "M", "mother", "daughter")
  )

# Load mother information and join
mother_info <- 
  read_csv("data/isala_pairs_mother_information.csv")

isala <- 
  isala %>%
  add_metadata(mother_info)

# Calculate Bray-Curtis dissimilarity index between mothers and their real daughters
betas <- 
  betas(isala, unique = FALSE) %>%
  filter(
    sample_type_pt_1 == "mother", 
    family_number_1 == family_number_2, beta != 0) %>%
  select(sample_1, code_1, sample_2, code_2, beta)

observed_beta_values <- betas %>% select(beta)

# Join beta info to isala
betas_info <- 
  betas %>%
  select(sample_2, beta) %>%
  rename(sample = sample_2)

isala <- isala %>% add_metadata(betas_info)

# Save updated tt object with beta info
#write_tidytacos(isala, "data/tt_isala_pairs_updated_betas")

# Create the permutation matrix 
num_permutations <- 10000
permuted_matrix <- matrix(nrow = 23, ncol = num_permutations)

for (i in 1:num_permutations) {
  mothers <- 
    isala %>% 
    filter_samples(isala$samples$sample_type_pt == "mother")
  daughters <- 
    isala %>% 
    filter_samples(isala$samples$sample_type_pt == "daughter")
  
  # Randomize mother-daughter pairing
  daughters$samples$mother_pt2 <- 
    sample(daughters$samples$mother_pt, replace = FALSE)
  names(daughters$samples$mother_pt2) <- 
    daughters$samples$mother_pt
  
  # Merge and compute betas
  isala2 <- 
    merge_tidytacos(mothers, daughters, taxon_identifier = taxon_name)
  permutedbeta <- 
    betas(isala2, unique = FALSE) %>%
    filter(sample_1 == mother_pt2_2) %>%
    select(beta)
  
  permuted_matrix[, i] <- permutedbeta$beta
}

# Save permutation matrix
write.csv(permuted_matrix, "data/isala_permutation_matrix.csv", row.names = FALSE)

# Load permutation matrix for visualization 

fin_matrix <- "data/isala_permutation_matrix.csv"
permuted_matrix <- read.csv(fin_matrix)

num_permutations <- 10000

# Density plot of mean beta values
tibble(mean_beta = colMeans(permuted_matrix)) %>%
  ggplot(aes(x = mean_beta)) +
  geom_density() +
  geom_rug(
    data = tibble(mean_beta_real = mean(
      observed_beta_values[[1]])), aes(x = mean_beta_real)) +
  xlim(c(0, 1))

# Summary stats
mean_random <- sum(colMeans(permuted_matrix)) / num_permutations
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
  geom_boxplot(
    width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_segment(
    data = tibble(
      mean_beta_real = mean(obs_similarity[[1]]), group = "Real pairs"),
    aes(x = group, xend = group, y = mean_beta_real, yend = mean_beta_real),
    color = "#00719A", size = 1) +
  geom_point(
    data = tibble(
      mean_beta = obs_similarity[[1]], group = "Real pairs"),
    aes(x = group, y = mean_beta), 
    color = "#00719A", size = 1.2, alpha = 0.5) +
  geom_segment(
    data = tibble(
      mean_beta_real = mean(obs_similarity[[1]])),
    aes(x = 1, xend = 2, y = mean_beta_real, yend = mean_beta_real),
    color = "#00719A", size = 1) +
  geom_text(
    data = tibble(
      mean_beta_real = mean(obs_similarity[[1]]), group = "Real pairs"),
    aes(x = group, y = mean_beta_real, label = round(mean_beta_real, 3)),
    color = "black", size = 4, hjust = -0.3) +
  theme_bw() +
  xlab("") +
  ylab("Mean Bray-Curtis Similarity Index (10K permutations)") +
  scale_x_discrete(
    limits = c("Real pairs", "Random pairs"),
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
ggsave(
  paste0(dout, "isala_permutation_test.svg"), violin_plot_comparison,
       units = "cm", width = 12, height = 24, dpi = 600)

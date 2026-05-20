# This script merges 16S and shotgun datasets, applies taxonomy corrections
# and prepares data for visualization.

# Dependencies: R, tidyverse, tidytacos, svglite

library(tidyverse)
library(tidytacos)
library(svglite)

# Define paths of input and output files
din_isaladmgs <- "data/tt_isala_shotgun_krakendb"
din_isala16s <- "data/tt_isala_16s_hq"
fin_metadata   <- "data/metadata.csv"
fin_lactobacilli <- "data/summary_lactobacillus_subgenera.txt"
dout <- "results/isala/quality_control"

# Make output folder 
if (!dir.exists(dout)) dir.create(dout, recursive = TRUE)

# Load metadata file
metadata <- 
  read_csv(fin_metadata) %>%
  rename(sample = bioSampleId)

# Load 16S sequencing data
isala_16s <- 
  read_tidytacos(din_isala16s) %>%
  mutate_samples(source = "16S") %>%
  add_metadata(metadata) %>%
  aggregate_taxa(rank = "genus") %>%
  select_taxa(c(taxon_id, genus))

# Load profiles based on shotgun data
isala_dmgs <- 
  read_tidytacos(din_isaladmgs) %>%
  set_rank_names(c(
    "domain", "phylum", "class", "order", "family", "genus", "species")) %>%
  mutate_taxa(species = str_remove(species, "^s_")) %>%
  add_metadata(metadata) %>%
  filter_samples(high_quality_16S == "TRUE") %>%
  select_samples(sample, sample_id)

# List Lactobacillus subgenus outline
lacto_map <-
  read_csv(fin_lactobacilli) %>%
  select(species, subgenus) %>%
  rename(taxon_name = species,
         genus = subgenus)

# Manually reclassify taxonomy of the tt object derived from shotgun data
isala_dmgs_genus <- 
  isala_dmgs %>%
  mutate_taxa(
    genus=str_replace(
      genus, "g_", "")
  ) %>%
  # Lactobacillus group
  mutate_taxa(
    genus = case_when(
      species == "Lactobacillus crispatus" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus delbrueckii" ~ 
        "Lactobacillus delbrueckii group",
      species == "Lactobacillus helveticus" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus johnsonii" ~ 
        "Lactobacillus gasseri group",
      species == "Lactobacillus acidophilus" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus taiwanensis" ~ 
        "Lactobacillus gasseri group",
      species == "Lactobacillus gasseri" ~ 
        "Lactobacillus gasseri group",
      species == "Lactobacillus acetotolerans" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus iners" ~ 
        "Lactobacillus iners group",
      species == "Lactobacillus amylovorus" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus paragasseri" ~ 
        "Lactobacillus gasseri group",
      species == "Lactobacillus bombicola" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus jensenii_A" ~ 
        "Lactobacillus jensenii group",
      species == "Lactobacillus gallinarum" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus jensenii" ~ 
        "Lactobacillus jensenii group",
      species == "Lactobacillus kefiranofaciens" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus amylolyticus" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus apis" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus equicursoris" ~ 
        "Lactobacillus delbrueckii group",
      species == "Lactobacillus melliventris" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus sp003692965" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus helsingborgensis" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus kullabergensis" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus sp002911475" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus pasteurii" ~ 
        "Lactobacillus pasteurii group",
      species == "Lactobacillus intestinalis" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus ultunensis" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus hominis" ~ 
        "Lactobacillus gasseri group",
      species == "Lactobacillus gigeriorum" ~ 
        "Lactobacillus pasteurii group",
      species == "Lactobacillus psittaci" ~ 
        "Lactobacillus jensenii group",
      species == "Lactobacillus delbrueckii_A" ~ 
        "Lactobacillus delbrueckii group",
      species == "Lactobacillus kitasatonis" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus hamsteri" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus kalixensis" ~ 
        "Lactobacillus crispatus group",
      species == "Lactobacillus delbrueckii_B" ~ 
        "Lactobacillus delbrueckii group",
      species == "Lactobacillus kimbladii" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus panisapium" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus sp000760615" ~ 
        "Lactobacillus apis group",
      species == "Lactobacillus rodentium" ~ 
        "Lactobacillus gasseri group",
      species == "Lactobacillus rodentium" ~ 
        "Lactobacillus gasseri group",
      TRUE ~ genus
    )) 

isala_dmgs_genus <-
  isala_dmgs_genus %>%
  aggregate_taxa(rank = "genus") %>%
  filter_taxa(!is.na(genus)) %>%
  add_taxon_name() 

isala_dmgs_genus <- 
  isala_dmgs_genus %>%
  mutate_samples(
    source = "Shotgun",
    participant = sample) 

isala_dmgs_genus <-
  isala_dmgs_genus %>%
  add_metadata(metadata) 

# Combine tt objects from 16S and DMGS data
combined_16s_mgs <- 
  merge_tidytacos(
    isala_16s, isala_dmgs_genus, taxon_identifier = genus)

combined_16s_mgs <-
  combined_16s_mgs %>% 
  set_rank_names("genus") 

# Visualization 

colors <- c(
  "#BEA5D5", "#8C564B", "#D62728", "#FF9896", "#FF7F0E",
  "#FFBB78", "#1F77B4", "#AEC7E8", "#2CA02C", "#FAF3DD"
)

combined_16s_mgs %>%
  tacoplot_stack(x = source, order_by = "family_number", n = 10) +
  facet_wrap(~ code, nrow = 6) +
  theme_bw() +
  xlab("Techniques") + 
  ylab("Relative abundance (%)") + 
  guides(fill = guide_legend(title = "Sub(genus)")) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(limits = c("16S", "Shotgun")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "italic"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.spacing.x = unit(1.5, "pt"),
    panel.grid = element_blank()
  )

# Save figure
ggsave(
  paste0(dout, "/comparison_per_participant.svg"), units = "cm", width = 40, 
  height = 20
)

# Add information about quality and remove samples from F20 onwards
metadata2 <-
  metadata %>%
  mutate(
    high_quality_shotgun = as.integer(str_remove(family_number, "F")) < 20
  )

write_csv(metadata2, "data/metadata_quality_info.csv")

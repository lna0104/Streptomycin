library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(ggh4x)
library(legendry)

source("R/plotting_nt.R")


# Load data
rrs_data <- read_csv("./output/rrs_filtered_output.csv")
rpsL_data <- read_csv("./output/rpsL_filtered_output.csv")

# Explore data
# Unique genomes in each dataset
rrs_genomes <- unique(rrs_data$accession_numbers)
rpsL_genomes <- unique(rpsL_data$accession_numbers)

# Genomes only in rrs_data
only_in_rrs <- setdiff(rrs_genomes, rpsL_genomes)
cat("Number of genomes only in rrs_data:", length(only_in_rrs), "\n")

# Genomes only in rpsL_data
only_in_rpsL <- setdiff(rpsL_genomes, rrs_genomes)
cat("Number of genomes only in rpsL_data:", length(only_in_rpsL), "\n")

# Genomes in both datasets
in_both <- intersect(rrs_genomes, rpsL_genomes)
cat("Number of genomes in both datasets:", length(in_both), "\n")

# Process data
rrs_res_list <- rrs_data |>
  filter(mutation_category == "present") |>
  distinct(accession_numbers) |>
  pull(accession_numbers)

rpsL_res_list <- rpsL_data |>
  filter(mutation_category == "present") |>
  distinct(accession_numbers) |>
  pull(accession_numbers)


# Get all unique accession numbers from both datasets
all_accessions <- union(rrs_data$accession_numbers, rpsL_data$accession_numbers)

# Create accession-to-genus mapping (one genus per accession)
meta_data <- bind_rows(
  rrs_data |> select(accession_numbers, species, genus),
  rpsL_data |> select(accession_numbers, species, genus)
) |>
  distinct(accession_numbers, .keep_all = TRUE)

# # Combine all species 
# all_species <- bind_rows(rrs_gene_list, rpsL_gene_list) %>%
#   distinct(species, genus)

# Create combined dataframe with mutation status and genus
combined_df <- tibble(
  accession_number = all_accessions,
  rrs_mutation = if_else(all_accessions %in% rrs_res_list, "present", "not_present"),
  rpsL_mutation = if_else(all_accessions %in% rpsL_res_list, "present", "not_present")
) |>
  left_join(meta_data, by = c("accession_number" = "accession_numbers"))


# Classify category
combined_species <- combined_df %>%
  mutate(
    category = case_when(
      rrs_mutation == "present" & rpsL_mutation == "present" ~ "both",
      rrs_mutation == "present" ~ "rrs",
      rpsL_mutation == "present" ~ "rpsL",
      TRUE ~ "none"
    )
  )

# Load gtdb taxonomy file from GTDB metadata
rrs_bacterial_taxonomy <- read_csv("./data/rrs_NCBI_taxonomy.csv", show_col_types = FALSE)
rpsL_bacterial_taxonomy <- read_csv("./data/NCBI_taxonomy.csv", show_col_types = FALSE) 

bacterial_taxonomy <- bind_rows(rrs_bacterial_taxonomy, rpsL_bacterial_taxonomy) |>
  distinct()

# Merge data with taxonomy
processed_data <- combined_species |>
  select(species, genus, category) |>
  left_join(bacterial_taxonomy, by = join_by(genus == genus))|>
  filter(!is.na(class)) |>
  group_by(class, phylum, category) |>
  summarise(
    n = n(),
    .groups = "drop"
  ) 

plot_class_counts_and_percents(processed_data, file_name = "./plots/rrs_rpsL_classes_genera.pdf")




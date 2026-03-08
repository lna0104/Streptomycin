t <- read_tsv("./data/bac120_metadata.tsv", show_col_types = FALSE) |>
  select(gtdb_taxonomy, ncbi_taxonomy)

library(dplyr)
library(tidyr)
library(stringr)

t1 <- t |>
  separate(
    gtdb_taxonomy,
    into = c("gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order", "gtdb_family", "gtdb_genus", "gtdb_species"),
    sep = ";\\s*",
    remove = TRUE,
    fill = "right",
    extra = "drop"
  ) |>
  separate(
    ncbi_taxonomy,
    into = c("ncbi_domain", "ncbi_phylum", "ncbi_class", "ncbi_order", "ncbi_family", "ncbi_genus", "ncbi_species"),
    sep = ";\\s*",
    remove = TRUE,
    fill = "right",
    extra = "drop"
  ) |>
  mutate(
    across(starts_with("gtdb_"), ~ str_remove(.x, "^[a-z]__")),
    across(starts_with("ncbi_"), ~ str_remove(.x, "^[a-z]__"))
  ) |>
  select(gtdb_genus, ncbi_genus) |>
  distinct()

ncbi_gtdb_species_map <- meta_data |>
  distinct(ncbi_species, species) |>
  filter(!is.na(ncbi_species), !is.na(species)) |>
  group_by(ncbi_species) |>
  filter(n_distinct(species) == 1) |> # remove one NCBI to multiple GTDB matching
  ungroup() |>
  group_by(species) |>
  filter(n_distinct(ncbi_species) == 1) |> # remove one GTDB to multiple NCBI matching
  ungroup()

t2 <- t1 |>
  group_by(gtdb_genus) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 1)


t3 <- t1 |>
  group_by(ncbi_genus) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 1)


ncbi_gtdb_species_map <- meta_data |>
  distinct(ncbi_species, species) |>
  filter(!is.na(ncbi_species), !is.na(species)) |>
  group_by(ncbi_species) |>
  filter(n_distinct(species) == 1) |> # remove one NCBI to multiple GTDB matching
  ungroup() |>
  group_by(species) |>
  filter(n_distinct(ncbi_species) == 1) |> # remove one GTDB to multiple NCBI matching
  ungroup()

tmp <- filtered_output |>
  left_join(bacterial_taxonomy, by = join_by(genus == genus)) |>
  distinct()

t <- read_csv("./results/rpsL_predicted_resistance.csv")
t1 <- t %>%  filter(mutations=="88R") %>% pull(species)
t2 <- t %>%  filter(mutations=="43R") %>% pull(species)
t3 <- t %>%  filter(mutations !="43R" & mutations != "88R")


t <- species_output |> filter(resistance == "resistant")
t1 <- filtered_output_gtdb |> filter(mutation_category == "present") |>
  pull(species) |>
  unique()

t2 <- output |>
  # step 1: merge multiple gene copies into one:
  group_by(species, genus, accession_numbers, mutation_name) |>
  summarise(
    mutation_category = ifelse(any(mutation_category == "present"),
                               "present",
                               ifelse(any(mutation_category == "possible"), "possible", "impossible")
    ),
    n_possible = sum(n_possible),
    .groups = "drop"
  ) |>
  # step 2: summarise across all screened mutations:
  group_by(species, genus, accession_numbers) |>
  summarise(
    resistance = any(mutation_category == "present"),
    evolvabilityI = sum(mutation_category == "possible"),
    evolvabilityII = sum(n_possible),
    .groups = "drop"
  ) 

length(t2 %>% filter(resistance == TRUE) %>% pull(species) %>% unique())
length(species_output %>% filter(resistance == "resistant") %>% pull(species) %>% unique())
length(species_data %>% filter(resistance == "resistant") %>% pull(species) %>% unique())
length(data_tree %>% filter(resistance == "resistant") %>% pull(label) %>% unique())
tmp <- species_output |>
  filter(resistance == "resistant") |>
  pull(species) |>
  setdiff(subtree$tip.label) %>% 
  as_tibble()

dtree <- as_tibble(subtree$tip.label)


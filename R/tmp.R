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

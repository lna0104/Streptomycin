library(readxl)
ncbi_gtdb_map <- read_excel("data/ncbi_vs_gtdb_bacteria.xlsx", sheet = 6)
ncbi_gtdb_map$gtdb_species <- sub(
  "^s__([^,]+?) [0-9.]+%.*",
  "\\1",
  ncbi_gtdb_map$`List of R226 species`
)



ncbi_gtdb_map$ncbi_species <- gsub(
  "\\[|\\]|^s__",
  "",
  ncbi_gtdb_map$`NCBI species`
)


write.csv(
  ncbi_gtdb_map |>
    select(ncbi_species, gtdb_species),
  "data/ncbi_gtdb_species.csv",
  row.names = FALSE
)

tmp <- ncbi_gtdb_species_map |>
  group_by(ncbi_species) |>
  summarise(n = n(), .groups = "drop")

tmp3 <- ncbi_gtdb_species_map |>
  group_by(species) |>
  summarise(n = n(), .groups = "drop")

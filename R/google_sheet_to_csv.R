# This script downloads data from google sheets, removes our own ALJElab mutations and some other mutations, 
# and saves the resulting tables as "reported_mutations.csv" and "report_references.csv".
# This file should not be included in the final version of the paper.

library(googlesheets4)
gs4_deauth() #function to inactive Google authorization
muts <- read_sheet("https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc/edit?usp=sharing",
                   sheet = "Data",
                   col_types = "c") |>
  filter(substr(Ref_code, 1, 4) != "ALJE") |>
  filter(Gene == "rpoB") |>
  filter(Origin %in% c("Isolate", "Lab mutant", "Construct")) |>
  mutate(Origin = ifelse(Origin %in% c("Lab mutant", "Construct"), "Lab-generated", Origin)) |>
  select(-Entry_by) |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)
write_csv(muts, "./data/reported_mutations.csv")


refs <- read_sheet("https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc/edit?usp=sharing",
                   sheet = "References",
                   col_types = "c") |>
  select(-Comments) |>
  mutate(DOI = sub("\\.$", "", Link)) |>
  select(-Link) |>
  arrange(Ref_code)

if (length(refs$Ref_code) != (length(unique(refs$Ref_code)))) {
  warning("Duplicated Ref_code entries detected!")
}

write_csv(refs, "./data/reports_references.csv")

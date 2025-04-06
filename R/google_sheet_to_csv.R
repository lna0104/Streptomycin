# This script downloads data from google sheets, removes our own ALJElab mutations and some other mutations, 
# and saves the resulting tables as "reported_mutations.csv" and "report_references.csv".
# This file should not be included in the final version of the paper.

library(googlesheets4)
library(tidyverse)

gs4_deauth() #function to inactive Google authorization
muts <- read_sheet("https://docs.google.com/spreadsheets/d/15M9eY1o_3uzA2zeR78SksPYUtUOOsT_e6XkRxPg6hA8/edit?gid=1074478647#gid=1074478647",
                   sheet = "Data",
                   col_types = "c") |>
  filter(substr(Ref_code, 1, 4) != "ALJE") |>
  filter(Gene == "rpsL" | Gene == "rrs") |>
  filter(Origin %in% c("Isolate", "Lab mutant", "Construct")) |>
  mutate(Origin = ifelse(Origin %in% c("Lab mutant", "Construct"), "Lab-generated", Origin)) |>
  select(-Entry_by) |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)
write_csv(muts, "./data/reported_mutations.csv")


refs <- read_sheet("https://docs.google.com/spreadsheets/d/15M9eY1o_3uzA2zeR78SksPYUtUOOsT_e6XkRxPg6hA8/edit?gid=992930929#gid=992930929",
                   sheet = "References",
                   col_types = "c") |>
  filter(Data_added=="yes")  |>
  select(-Comments) |>
  mutate(DOI = sub("\\.$", "", Link)) |>
  select(-Link) |>
  arrange(Ref_code)

filtered_refs <- refs |> filter(Ref_code %in% muts$Ref_code)

if (length(filtered_refs$Ref_code) != (length(unique(filtered_refs$Ref_code)))) {
  warning("Duplicated Ref_code entries detected!")
}


write_csv(refs, "./data/reports_references.csv")

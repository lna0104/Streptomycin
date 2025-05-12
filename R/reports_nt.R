#' Summarise reported mutations
#'
#' @param muts a data frame providing info on mutations reported in the literature
#' @param min_n_species a numeric vector showing the minimum number of species which a mutation has been reported in
#' @param file_name a path that the summary should be saved in
#'
#' @return null
#' @export
#'
#' @examples summarise_reported_mutation (muts, min_n_species, "./output/myfilename.txt")
summarise_reported_mutations_nt <- function(muts, file_name, subtitle = "") {
  
  # studies:
  n_studies <- muts |>
    pull(Ref_code) |>
    unique() |>
    length()
  
  # positions with mutations:
  n_positions <- muts %>%
    filter(!is.na(Nt_pos_Ecoli), !is.na(Nt_mutation)) %>%
    distinct(Nt_pos_Ecoli) %>% 
    pull() |>
    length()
  
  # Nt mutations: 
  n_Nt_mutations <- muts %>%
    filter(!is.na(Nt_pos_Ecoli), !is.na(Nt_mutation)) %>%
    mutate(mutation_name = paste(Nt_pos_Ecoli, Nt_mutation, sep = "_")) %>% 
    distinct(mutation_name) %>% 
    pull() |>
    length()
  
  #reported species
  n_species <- muts |>
    filter(!is.na(Nt_pos_Ecoli), !is.na(Nt_mutation)) |>
    pull(Species) |> 
    unique() |>
    length()
  
  report <- paste0(
    "Summary of reported mutations ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "Number of studies: ", n_studies, "\n",
    "Number of species: ", n_species, "\n",
    "Total number of mutations: ", nrow(muts), "\n",
    "Nucleotide positions with mutations: ", n_positions, "\n",
    "Number of unique nucleotide mutations: ", n_Nt_mutations, "\n"
  )
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}

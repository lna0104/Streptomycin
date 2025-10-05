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

#' Summarise mutation screen results
#'
#' @param filtered_output filtered output from mutation screen
#' @param file_name a path that the summary should be saved in
#' @param subtitle Additional information to add to the title
#' @param target_gene name of target gene
#' #'
#' @return null
#' @export
#'
#' @examples summarise_mutation_screen(final_output, filtered_output, "./output/myfilename.txt")
summarise_mutation_screen_nt <- function(filtered_output, target_gene, file_name, subtitle = ""){
   
  # species-level summary of filtered output:
  species_output <- get_species_output_nt(filtered_output)
  
  # mutations screened:
  n_screened_mutations <- filtered_output |>
    pull(mutation_name) |>
    unique() |> 
    length()
  n_sites <- filtered_output |>
    pull(Nt_pos_Ecoli) |>
    unique() |>
    length()
  
  # prediction of intrinsic resistance:
  
  #number of mutations present:
  n_present_mutations <- filtered_output |>
    filter(mutation_category == "present") |>
    nrow()
  
  # number of unique mutations present:
  n_unique_mutations <- filtered_output |>
    filter(mutation_category == "present") |>
    pull(mutation_name) |>
    unique() |>
    length()
  
  # number of species predicted to be resistant:
  n_resistant_species <- species_output |>
    filter(resistance == "resistant") |>
    nrow()
  
  # fraction of species predicted to be resistant:
  n_species <- filtered_output |> pull(species) |> unique() |> length()
  f_resistant_species <- n_resistant_species / n_species
  
  # number of species with more than one resistance mutation:
  n_multiresistant_species <- filtered_output |>
    filter(mutation_category == "present") |>
    group_by(species) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n > 1) |>
    nrow()
  
  # determine multicopy species:
  filtered_output_seqs <- filtered_output |>
    select(accession_numbers, species, gene_copy, target_length, alig_score, core_dist) |>
    distinct()
  multicopy_species <- filtered_output_seqs |>
    group_by(accession_numbers) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n > 1L) |>
    pull(accession_numbers)
  
  # number of species where one gene copy confers resistance and one doesn't:
  n_hetero_resistance <- filtered_output |>
    filter(accession_numbers %in% multicopy_species) |>
    group_by(species, accession_numbers, gene_copy) |>
    summarise(resistant_copy = any(mutation_category == "present"), .groups = "drop") |>
    group_by(species, accession_numbers) |>
    summarise(one_but_not_all_resistant = (any(resistant_copy) & (!all(resistant_copy))), .groups = "drop") |>
    filter(one_but_not_all_resistant) |>
    nrow()
  

  report <- paste0(
    "Summary of the mutation screen ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "Screened mutations:\n",
    "    Number of distinct mutations screened: ", n_screened_mutations, "\n",
    "    Number of nucleotide sites: ", n_sites, "\n",
    "Predicted intrinsic resistance:\n",
    "    Total number of resistance mutations present across species: ", n_present_mutations, "\n",
    "    Number of unique mutations present across species: ", n_unique_mutations, "\n",
    "    Number of resistant species: ", n_resistant_species, "\n",
    "    Percentage of resistant species: ", format(round(100 * f_resistant_species, 2), nsmall = 2), "%\n",
    "    Number of species with multiple resistance mutations: ", n_multiresistant_species, "\n",
    "    Number of species where one ", target_gene, " copy confers resistance and one does not: ", n_hetero_resistance, "\n"
  )
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}


summarise_phylogenetics_nt <- function(subtree, species_output, sample_n = NULL, file_name, subtitle = "") {
  
  phylo_signals <- get_phylosignals_nt(subtree, species_output, sample_n = sample_n)
  subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))
  
  n_resistant <- species_output |>
    filter(resistance == "resistant") |>
    pull(species) |>
    base::intersect(subtree$tip.label) |>
    length()
  
  if (is.null(sample_n))
    sample_n <- "all"
  
  report <- paste0(
    "Summary of the phylogenetic analyses ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "Phylogenetic tree:\n",
    "    Number of species: ", Ntip.phylo(subtree), "\n",
    "    Number of species predicted to be resistant in tree: ", n_resistant, "\n",
    "Number of species sampled from tree when calculating phylogenetic signals: ", sample_n, "\n",
    "Phylogenetic signal in predicted resistance:\n",
    "    Test: permutation test of mean phylogenetic distance of resistant species\n",
    "    Number of permutations: ", length(phylo_signals$permtest_resistance$mean_distance_tips_permutated), "\n",
    "    p-value: ", phylo_signals$permtest_resistance$P, "\n")
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}


summarise_conservation_nt <- function(cons, target_gene, file_name, subtitle = "") {
  
  report <- paste0(
    "Summary of the nucleotide conservation analyses ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "Mean Hamming distance to E. coli across all sequences:\n",
    "    Number of sequences: ", nrow(cons$hamming_Ecoli), "\n",
    "    Mean across nucleotide positions: ", mean(cons$means$hamming_Ecoli, na.rm = TRUE),  "\n",
    "    Max across nucleotide positions: ", max(cons$means$hamming_Ecoli, na.rm = TRUE),  "\n",
    "    Min across nucleotide positions: ", min(cons$means$hamming_Ecoli, na.rm = TRUE),  "\n",
    "Mean Hamming distance across randomly sampled pairs of", target_gene, " sequences:\n",
    "    Number of sequence pairs: ", nrow(cons$hamming_rnd), "\n",
    "    Mean across nucleotide positions: ", mean(cons$means$hamming_rnd, na.rm = TRUE),  "\n",
    "    Max across nucleotide positions: ", max(cons$means$hamming_rnd, na.rm = TRUE),  "\n",
    "    Min across nucleotide positions: ", min(cons$means$hamming_rnd, na.rm = TRUE),  "\n")
    # "Mean Grantham distance to E. coli across all sequences:\n",
    # "    Number of sequences: ", nrow(cons$grantham_Ecoli), "\n",
    # "    Mean across AA positions: ", mean(cons$means$grantham_Ecoli, na.rm = TRUE),  "\n",
    # "    Max across AA positions: ", max(cons$means$grantham_Ecoli, na.rm = TRUE),  "\n",
    # "    Min across AA positions: ", min(cons$means$grantham_Ecoli, na.rm = TRUE),  "\n",
    # "Mean Grantham distance across randomly sampled pairs of", target_gene, " sequences:\n",
    # "    Number of sequence pairs: ", nrow(cons$grantham_rnd), "\n",
    # "    Mean across AA positions: ", mean(cons$means$grantham_rnd, na.rm = TRUE),  "\n",
    # "    Max across AA positions: ", max(cons$means$grantham_rnd, na.rm = TRUE),  "\n",
    # "    Min across AA positions: ", min(cons$means$grantham_rnd, na.rm = TRUE),  "\n")
  
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}



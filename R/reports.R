#' Summarise reported mutations
#'
#' @param muts a data frame providing info on rpoB mutations reported in the literature
#' @param min_n_species a numeric vector showing the minimum number of species which a mutation has been reported in
#' @param file_name a path that the summary should be saved in
#'
#' @return null
#' @export
#'
#' @examples summarise_reported_mutation (muts, min_n_species, "./output/myfilename.txt")
summarise_reported_mutations <- function(muts, file_name, subtitle = "") {
  
  # stats before filtering:
  # studies:
  n_studies <- muts |>
    pull(Ref_code) |>
    unique() |>
    length()
  
  # positions with mutations:
  n_positions <- muts %>%
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) %>%
    distinct(AA_pos_Ecoli) %>% 
    pull() |>
    length()
  
  # AA mutations: 
  n_AA_mutations <- muts %>%
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) %>%
    mutate(mutation_name = paste(AA_pos_Ecoli, AA_mutation, sep = "_")) %>% 
    distinct(mutation_name) %>% 
    pull() |>
    length()
  
  #reported species
  n_species <- muts |>
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) |>
    pull(Species) |> 
    unique() |>
    length()
  
  # stats after filtering:
  mutation_list_reports <- filter_mutations(muts,
                                            min_n_species = globsets$min_n_species, 
                                            min_n_studies = globsets$min_n_studies)
  muts_F <- muts |>
    semi_join(mutation_list_reports, by = join_by(AA_pos_Ecoli, AA_mutation))
  # studies:
  n_studies_F <- muts_F |>
    pull(Ref_code) |>
    unique() |>
    length()
  
  # positions with mutations:
  n_positions_F <- muts_F %>%
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) %>%
    distinct(AA_pos_Ecoli) %>% 
    pull() |>
    length()
  
  # AA mutations: 
  n_AA_mutations_F <- muts_F %>%
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) %>%
    mutate(mutation_name = paste(AA_pos_Ecoli, AA_mutation, sep = "_")) %>% 
    distinct(mutation_name) %>% 
    pull() |>
    length()
  
  #reported species
  n_species_F <- muts_F |>
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) |>
    pull(Species) |> 
    unique() |>
    length()
  
  # Results relating to MTB:
  # number of filtered mutations reported in clinical MTB isolates:
  mtb_clinical <- muts |>
    filter(Species == "Mycobacterium tuberculosis" & Origin == "Isolate") |>
    distinct(AA_pos_Ecoli, AA_mutation)
  
  n_clinicalMTB_and_F <- mutation_list_reports |>
    semi_join(mtb_clinical, by = join_by(AA_pos_Ecoli, AA_mutation)) |>
    nrow()
  
  n_F_not_clinicalMTB <- mutation_list_reports |>
    anti_join(mtb_clinical, by = join_by(AA_pos_Ecoli, AA_mutation)) |>
    nrow()
  
  mutation_list_reports_no_MTB <- muts |>
    filter(Species != "Mycobacterium tuberculosis") |>
    filter_mutations(min_n_species = globsets$min_n_species, 
                     min_n_studies = globsets$min_n_studies)
  
  n_filtered_mutations_no_MTB <- nrow(mutation_list_reports_no_MTB)
  n_MTB_predicted_without_MTB <- mutation_list_reports_no_MTB |>
    semi_join(mtb_clinical) |>
    nrow()
  
  report <- paste0(
    "Summary of reported mutations ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "All mutations:\n",
    "    Number of studies: ", n_studies, "\n",
    "    Number of species: ", n_species, "\n",
    "    Total number of mutations: ", nrow(muts), "\n",
    "    Amino acid positions with mutations: ", n_positions, "\n",
    "    Number of unique amino acid mutations: ", n_AA_mutations, "\n",
    "Filtered mutations:\n",
    "    Number of studies: ", n_studies_F, "\n",
    "    Number of species: ", n_species_F, "\n",
    "    Total number of mutations: ", nrow(muts_F), "\n",
    "    Amino acid positions with mutations: ", n_positions_F, "\n",
    "    Number of unique amino acid mutations: ", n_AA_mutations_F, "\n",
    "Mutations in Mycobacterium tuberculosis (MTB) clinical isolates:\n",
    "    Number of filtered mutations also reported in MTB: ", n_clinicalMTB_and_F, "\n",
    "    Number of filtered mutations collated without MTB data: ", n_filtered_mutations_no_MTB, "\n",
    "    Number of filtered mutations collated without MTB data also reported in MTB: ", n_MTB_predicted_without_MTB, "\n"
  )
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}


#' Summarise target sequences screened
#' 
#' @param genome_summary summaries of downloaded genomes
#' @param final_output raw output from mutation screen
#' @param filtered_output filtered output from mutation screen
#' @param min_seq_length minimum sequence length
#' @param min_alig_score minimum alignment score of gene sequences
#' @param max_core_dist maximum Levenshtein distince of gene sequences to E. coli
#' @param file_name a path that the summary should be saved in
#' @param target_gene name of target gene
#' @param subtitle Additional information to add to the title
#' 
#' @return null
#' @export
#'
#' @examples summarise_target_sequences (muts, min_n_species, "./output/myfilename.txt")
summarise_target_sequences <- function(genome_summary, 
                                       final_output, 
                                       filtered_output, 
                                       min_seq_length,
                                       min_alig_score,
                                       max_core_dist,
                                       file_name, 
                                       target_gene,
                                       subtitle = "") {
  # stats before filtering:
  
  #number of genomes attempted to download:
  n_genomes_summary <- genome_summary %>% length()
  
  final_output_seqs <- final_output |>
    select(accession_numbers, species, gene_copy, target_length, alig_score, core_dist) |>
    distinct()
  #number of genomes with extracted gene sequences:
  n_analysed_genomes <- final_output_seqs$accession_numbers %>% unique %>% length()
  # number of genomes where no gene sequence could be extracted:
  n_unsuccessful_genomes <- n_genomes_summary - n_analysed_genomes
  # number of species among downloaded genomes:
  n_species <- final_output_seqs$species |> unique() |> length()
  # number of gene sequences:
  n_gene_sequences <- nrow(final_output_seqs)
  # number of species with >1 gene copy:
  copies_per_sequence <- final_output_seqs |>
    group_by(accession_numbers, species) |>
    summarise(n = n(), .groups = "drop")
  n_multicopy_species <-  copies_per_sequence |>
    filter(n > 1L) |> 
    nrow()
  # maximum number of gene sequences per species:
  n_max_gene_copy <- max(copies_per_sequence$n)
  
  # lengths:
  min_length <- min(final_output_seqs$target_length)
  median_length <- median(final_output_seqs$target_length)
  mean_length <- mean(final_output_seqs$target_length)
  max_length <- max(final_output_seqs$target_length)
  
  #aligning scores and dist
  min_alig_score <- min(final_output$alig_score)
  max_alig_score <- max(final_output$alig_score)
  max_dist_from_ref <- max(final_output$core_dist)
  min_dist_from_ref <- min(final_output$core_dist)
  
  # stats after filtering:
  
  filtered_output_seqs <- filtered_output |>
    select(accession_numbers, species, gene_copy, target_length, alig_score, core_dist) |>
    distinct()
  #number of genomes with extracted gene sequences:
  n_analysed_genomes_F <- filtered_output_seqs$accession_numbers %>% unique %>% length()
  # number of species among downloaded genomes:
  n_species_F <- filtered_output_seqs$species |> unique() |> length()
  # number of gene sequences:
  n_gene_sequences_F <- nrow(filtered_output_seqs)
  # number of species with >1 gene copy:
  copies_per_sequence_F <- filtered_output_seqs |>
    group_by(accession_numbers, species) |>
    summarise(n = n(), .groups = "drop")
  n_multicopy_species_F <-  copies_per_sequence_F |>
    filter(n > 1L) |> 
    nrow()
  # maximum number of gene sequences per species:
  n_max_gene_copy_F <- max(copies_per_sequence_F$n)
  
  report <- paste0(
    "Summary of screened target ", target_gene, " sequences ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "Target sequence statistics before filtering:\n",
    "    Number of genomes in search results: ", n_genomes_summary, "\n",
    "    Number of genomes with extracted ", target_gene, " sequence(s): ", n_analysed_genomes, "\n",    "    Number of genomes where no ", target_gene ," sequence could be extracted: ", n_unsuccessful_genomes, "\n",
    "    Number of species with downloaded genomes: ", n_species, "\n",
    "    Total number of ", target_gene, " sequences: ", n_gene_sequences, "\n",
    "    Number of genomes with more than one ", target_gene, " sequence: ", n_multicopy_species, "\n",
    "    Maximum number of ", target_gene, " sequences per genome: ", n_max_gene_copy, "\n",
    "Filtering statistics:\n",
    "    Minimum sequence length: ", min_length, "\n",
    "    Median sequence length: ", format(median_length, nsmall = 2), "\n",
    "    Mean sequence length: ", format(mean_length, nsmall = 2), "\n",
    "    Maximum sequence length: ", max_length, "\n",
    "    Minimum aligning score: ", format(min_alig_score, nsmall = 2), "\n",
    "    Maximum aligning score: ", format(max_alig_score, nsmall = 2), "\n",
    "    Minimum distance from reference ", target_gene, " core: ", min_dist_from_ref, "\n",
    "    Maximum distance from reference ", target_gene, " core: ", max_dist_from_ref, "\n",
    "Filters applied: ", "\n",
    "    Minimum sequence length: ", min_seq_length, "\n",
    "    Minimum alignment score to E. coli ", target_gene, ": ", min_alig_score, "\n",
    "    Maximum core distance to E. coli ", target_gene, ": ", max_core_dist, "\n",
    "Target sequence statistics after filtering:\n",
    "    Number of genomes with extracted ", target_gene, " sequence(s): ", n_analysed_genomes_F, "\n",
    "    Number of species with downloaded genomes: ", n_species_F, "\n",
    "    Total number of ", target_gene, " sequences: ", n_gene_sequences_F, "\n",
    "    Number of genomes with more than one ", target_gene, " sequence: ", n_multicopy_species_F, "\n",
    "    Maximum number of ", target_gene, " sequences per genome: ", n_max_gene_copy_F, "\n"
  )
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
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
summarise_mutation_screen <- function(filtered_output, target_gene, file_name, subtitle = ""){
   
  # species-level summary of filtered output:
  species_output <- get_species_output(filtered_output)
  
  # mutations screened:
  n_screened_mutations <- filtered_output |>
    pull(mutation_name) |>
    unique() |> 
    length()
  n_sites <- filtered_output |>
    pull(AA_pos_Ecoli) |>
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
  
  # quantiles for evolvability:
  quant_evolvabilityI <- species_output |>
    pull(evolvabilityI) |>
    quantile(p = c(0, 0.025, 0.5, 0.975, 1))
  
  quant_evolvabilityII <- species_output |>
    pull(evolvabilityII) |>
    quantile(p = c(0, 0.025, 0.5, 0.975, 1))
  
  # theoretical evolvabilities:
  theoretical_evolvabilities <- get_theoretical_evolvabilities(filtered_output)

  # associations between intrinsic resistance and evolvabilities:
  corr_evolvabilityI_vs_II <- cor(species_output$evolvabilityI, species_output$evolvabilityII)
  t_test_evolvabilityI <- t.test(species_output |> filter(resistance == "resistant") |> pull(evolvabilityI),
                                 species_output |> filter(resistance == "susceptible") |> pull(evolvabilityI))
  t_test_evolvabilityII <- t.test(species_output |> filter(resistance == "resistant") |> pull(evolvabilityII),
                                  species_output |> filter(resistance == "susceptible") |> pull(evolvabilityII))
  
  report <- paste0(
    "Summary of the mutation screen ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "Screened mutations:\n",
    "    Number of distinct mutations screened: ", n_screened_mutations, "\n",
    "    Number of amino acid sites: ", n_sites, "\n",
    "Predicted intrinsic resistance:\n",
    "    Total number of resistance mutations present across species: ", n_present_mutations, "\n",
    "    Number of unique mutations present across species: ", n_unique_mutations, "\n",
    "    Number of resistant species: ", n_resistant_species, "\n",
    "    Percentage of resistant species: ", format(round(100 * f_resistant_species, 2), nsmall = 2), "%\n",
    "    Number of species with multiple resistance mutations: ", n_multiresistant_species, "\n",
    "    Number of species where one ", target_gene, " copy confers resistance and one does not: ", n_hetero_resistance, "\n",
    "Evolvability I (number of AA mutations that a species can mutate to):\n",
    "    Range: ", quant_evolvabilityI[1], "...", quant_evolvabilityI[5], "\n",
    "    95% inter-quantile range: ", quant_evolvabilityI[2], "...", quant_evolvabilityI[4], "\n",
    "    Median: ", quant_evolvabilityI[3], "\n",
    "    Theoretical range: ", theoretical_evolvabilities$min_evolvabilityI, "...",
                               theoretical_evolvabilities$max_evolvabilityI, "\n",
    "Evolvability II (number of mutations that can produce a resistance mutation):\n",
    "    Range: ", quant_evolvabilityII[1], "...", quant_evolvabilityII[5], "\n",
    "    95% inter-quantile range: ", quant_evolvabilityII[2], "...", quant_evolvabilityII[4], "\n",
    "    Median: ", quant_evolvabilityII[3], "\n",
    "    Theoretical range: ", theoretical_evolvabilities$min_evolvabilityII, "...",
                               theoretical_evolvabilities$max_evolvabilityII, "\n",
    "Associations between evolvabilities (phylogenetically uncontrolled):\n",
    "    Evolvability I vs. evolvability II: r=", format(corr_evolvabilityI_vs_II, digits = 3), "\n",
    "    Intrinsic resistance vs. evolvability I: p=", format(t_test_evolvabilityI$p.value, digits = 3), 
      " (", t_test_evolvabilityI$method, ", t=", format(t_test_evolvabilityI$statistic, digits = 3), ")\n",
    "    Intrinsic resistance vs. evolvability II: p=", format(t_test_evolvabilityII$p.value, digits = 3), 
    " (", t_test_evolvabilityII$method, ", t=", format(t_test_evolvabilityII$statistic, digits = 3), ")\n"
  )
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}


summarise_phylogenetics <- function(subtree, species_output, sample_n = NULL, file_name, subtitle = "") {
  
  phylo_signals <- get_phylosignals(subtree, species_output, sample_n = sample_n)
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
    "    p-value: ", phylo_signals$permtest_resistance$P, "\n",
    "Phylogenetic signal in evolvability I (number of evolvable AA mutations):\n",
    "    Pagel's lambda: ", format(phylo_signals$lambda_evolvabilityI$lambda, digits = 3), "\n",
    "    p(lambda): ", format(phylo_signals$lambda_evolvabilityI$P, digits = 3), "\n",
    "    Blomberg's K: ", format(phylo_signals$K_evolvabilityI$K, digits = 3), "\n",
    "    p(K): ", format(phylo_signals$K_evolvabilityI$P, digits = 3), "\n",
    "Phylogenetic signal in evolvability II (number of nt mutations producing AA resistance mutations):\n",
    "    Pagel's lambda: ", format(phylo_signals$lambda_evolvabilityII$lambda, digits = 3), "\n",
    "    p(lambda): ", format(phylo_signals$lambda_evolvabilityII$P, digits = 3), "\n",
    "    Blomberg's K: ", format(phylo_signals$K_evolvabilityII$K, digits = 3), "\n",
    "    p(K): ", format(phylo_signals$K_evolvabilityII$P, digits = 3), "\n")
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}


summarise_conservation <- function(cons, target_gene, file_name, subtitle = "") {
  
  report <- paste0(
    "Summary of the amino acid conservation analyses ", subtitle, "\n",
    "Date: ", Sys.time(), "\n",
    "--------------------------------------------------------------------\n\n",
    "Mean Hamming distance to E. coli across all sequences:\n",
    "    Number of sequences: ", nrow(cons$hamming_Ecoli), "\n",
    "    Mean across AA positions: ", mean(cons$means$hamming_Ecoli, na.rm = TRUE),  "\n",
    "    Max across AA positions: ", max(cons$means$hamming_Ecoli, na.rm = TRUE),  "\n",
    "    Min across AA positions: ", min(cons$means$hamming_Ecoli, na.rm = TRUE),  "\n",
    "Mean Hamming distance across randomly sampled pairs of", target_gene, " sequences:\n",
    "    Number of sequence pairs: ", nrow(cons$hamming_rnd), "\n",
    "    Mean across AA positions: ", mean(cons$means$hamming_rnd, na.rm = TRUE),  "\n",
    "    Max across AA positions: ", max(cons$means$hamming_rnd, na.rm = TRUE),  "\n",
    "    Min across AA positions: ", min(cons$means$hamming_rnd, na.rm = TRUE),  "\n",
    "Mean Grantham distance to E. coli across all sequences:\n",
    "    Number of sequences: ", nrow(cons$grantham_Ecoli), "\n",
    "    Mean across AA positions: ", mean(cons$means$grantham_Ecoli, na.rm = TRUE),  "\n",
    "    Max across AA positions: ", max(cons$means$grantham_Ecoli, na.rm = TRUE),  "\n",
    "    Min across AA positions: ", min(cons$means$grantham_Ecoli, na.rm = TRUE),  "\n",
    "Mean Grantham distance across randomly sampled pairs of", target_gene, " sequences:\n",
    "    Number of sequence pairs: ", nrow(cons$grantham_rnd), "\n",
    "    Mean across AA positions: ", mean(cons$means$grantham_rnd, na.rm = TRUE),  "\n",
    "    Max across AA positions: ", max(cons$means$grantham_rnd, na.rm = TRUE),  "\n",
    "    Min across AA positions: ", min(cons$means$grantham_rnd, na.rm = TRUE),  "\n")
  
  write_file(report, file_name)
  cat(report)
  cat(paste0("\n\nThis summary has been saved in file ", file_name, ".\n"))
  return(invisible(NULL))
}



render_summary <- function(summaries, preamble, summaries_path, show_all_dates = FALSE) {
  
  lines <- vector(mode = "list", length = length(summaries) + 2)
  lines[[1]] <- readLines(preamble)
  
  for(i in 2:(length(summaries) + 1)) {
    lines[[i]] <- readLines(paste0(summaries_path, "/summary_", summaries[i - 1], ".txt"))
    lines[[i]][1] <- paste("##", lines[[i]][1])
    if (show_all_dates) {
      lines[[i]][2] <- paste0("*", lines[[i]][2], "*")
      
    } else {
      lines[[i]][2] <- ""
    }
    lines[[i]] <- lines[[i]][-3]
    lines[[i]][3:length(lines[[i]])] <- sub("^(\\S)", "* \\1", lines[[i]][3:length(lines[[i]])])
    lines[[i]][3:length(lines[[i]])] <- sub("^ {4}", "    + ", lines[[i]][3:length(lines[[i]])])
    lines[[i]] <- c("", lines[[i]][1], "", lines[[i]][2:length(lines[[i]])])
  }
  lines[[length(summaries) + 2]] <- c("", "## Session Info", "", pander::pander_return(sessionInfo()))
  lines <- unlist(lines)
  write_lines(lines, paste0(summaries_path, "/summary.qmd"))
  quarto_render(paste0(summaries_path, "/summary.qmd"))
}


make_table_intrinsic_resistance <- function(filtered_output, file_name) {
  filtered_output |>
    filter(mutation_category == "present") |>
    mutate(mutation_name = str_replace(mutation_name, "_", "")) |>
    group_by(accession_numbers, species) |>
    summarise(mutations = paste(mutation_name, collapse = ", "), .groups = "drop") |>
    rename(accession_number = accession_numbers) |>
    select(species, accession_number, mutations) |>
    arrange(species) |>
    write_csv(file_name)
}

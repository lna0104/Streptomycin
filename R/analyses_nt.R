#' Prepare a (possibly filtered) table of mutations to be screened:
#'
#' @param muts a data frame providing info on rrs mutations reported in the literature
#' @param min_n_species a numeric vector providing the minimum number of species which a mutation has been reported in
#' @param min_n_studies a numeric vector providing the minimum number studies which a mutation has been reported in
#' @param origin Origin of mutations to be included, e.g., "Isolate" or "Lab mutant". 
#' The default value of "ANY" means that all origins will be included. 
#'
#' @return a data frame providing info on filtered reported mutations in the literature(reliable mutations)
#' @export
#'
#' @examples filter_mutations_nt(muts, min_n_species = 5)
filter_mutations_nt <- function(muts, 
                                min_n_species = 1, 
                                min_n_studies = 1, 
                                origin = "ANY") {
  
  mutation_list <- muts |>
    filter(!is.na(Nt_pos_Ecoli), !is.na(Nt_original), !is.na(Nt_mutation)) 
  # filter(is.na(Warning) | Warning == "AA_pos inconsistent with AA_pos_Ecoli")
  
  if (origin == "ANY") {
    origin <- muts |> pull(Origin) |> unique()
  } 
  
  mutation_list_species <- mutation_list |>
    select(Species, Nt_pos_Ecoli, Nt_mutation) |>
    distinct() |>
    group_by(Nt_pos_Ecoli, Nt_mutation) |>
    summarise(n_species = n(), .groups = "drop")
  
  mutation_list_studies <- mutation_list |>
    select(Ref_code, Nt_pos_Ecoli, Nt_mutation) |>
    distinct() |>
    group_by(Nt_pos_Ecoli, Nt_mutation) |>
    summarise(n_studies = n(), .groups = "drop")
  
  mutation_list <- mutation_list_studies |>
    full_join(mutation_list_species, by = join_by(Nt_pos_Ecoli, Nt_mutation)) |>
    filter(n_studies >= min_n_studies) |>
    filter(n_species >= min_n_species) |>
    select(Nt_pos_Ecoli, Nt_mutation)
  
  return(mutation_list)
}

#' Screen a set of target sequences for existence and evolvability over list of mutations
#'
#' @param target_sequences A DNAStringSet object containing all target gene sequences
#' @param mutation_list A data frame containing all the mutations to be screened for
#'
#' @return A list of either data frames (if screen for a sequence was successful) or error messages (if not).
#'
#'
screen_target_sequences_nt <- function(target_sequences, reference_Ecoli, mutation_list, target_gene, n_workers) {
  # Ensure required directories exist
  if (!dir.exists("./output/alignments_nt")) dir.create("./output/alignments_nt", recursive = TRUE)
  if (!dir.exists("./output/progress_nt")) dir.create("./output/progress_nt", recursive = TRUE)
  
  plan(multisession, workers = n_workers) # to parallelise the function
  final_output <- future_lapply(1:length(target_sequences), function(i) {
    cat(paste0("Working on species #", i, "/", length(target_sequences), "...\n")) #notify when the analysis of each sequence starts
    output <- tryCatch(
      check_out_mutation_nt(mutation_list,              # compare target sequences with reference sequence to identify the nucleotide at each position
                            target_sequences[[i]], 
                            reference_Ecoli,
                            target_gene,
                            alig_file = paste0("./output/alignments_nt/", names(target_sequences)[[i]], ".RData")),
      error = function(cond) {
        message("  It seems there is a problem with this sequence. Here's the original error message:") # gives error but keeps the analysis moves on
        message(conditionMessage(cond)) #show what is the error
        return(list(error = conditionMessage(cond)))
      }
    )
    if(!is.null(output)) {
      output$accession_numbers <- extract_accession_number(names(target_sequences)[i], target_gene) #insert a column for each target gene assembly accession number
      output$species <- extract_species_name(names(target_sequences)[i], target_gene)  #insert a column for each species name
      output$target_name <- names(target_sequences)[i]
      output$gene_copy <- extract_gene_copy_number(names(target_sequences)[i], target_gene) #insert a column for each target gene copy number
      output$target_length <- length(target_sequences[[i]])
      if (is.data.frame(output)) {
        output <- select(output, accession_numbers:target_length, alig_score, core_dist, everything())
      }
    }
    
    saveRDS(i, file = paste0("./output/progress_nt/sequence_", i, "_out_of_", length(rrs_target_sequences), "_complete.rds"))
    return(output)
  })
  return(final_output)
}

#' Title: check out mutation nt
#'
#' @param mutations_list a data frame providing all nucleotide substitutions reported in Google sheet based on coordinations
#' @param target_sequences a list of multiple DNAString objects (gene sequences retrieved from different bacterial species) 
#' @param reference_Ecoli a DNAString object (gene sequence of Escherichia coli MG1655)
#' @param target_gene a name of target gene
#' @return A data frame of all nucleotide substitutions within all target sequences and check if there is a nucleotide substantiation and whether nucleotide substantiation is new or has been reported before or there is a different nucleotide in the original position 
#' @export
#'
#' @examples check_out_mutation(mutations_list, target_sequence, reference_Ecoli, target_gene)
check_out_mutation_nt <- function(mutations_list, 
                                  target_sequence, 
                                  reference_Ecoli,
                                  target_gene,
                                  alig_file = NULL) {
  if (is.null(target_sequence) || length(target_sequence) == 0) {
    stop(sprintf('Empty "%s" sequence.', target_gene))
  }
  coords_alig <- getCoordinatesNt(target_sequence, reference_Ecoli, aligOutput = TRUE)
  if (!is.null(alig_file)) {
    save(coords_alig, file = alig_file)
  }
  coords <- coords_alig$coordinates
  # protein_target <- Biostrings::translate(target_sequence)
  # protein_Ecoli <- Biostrings::translate(reference_Ecoli)
  
  result <- mutations_list
  result$Nt_target <- NA # what is the nucleotide at that position in focal species
  result$n_possible <- NA # number of single nt mutations in sequence that would produce the mutation
  for (j in 1:nrow(mutations_list)) {
    Nt_pos <- ALJEbinf::translateCoordinate(mutations_list$Nt_pos_Ecoli[j], 
                                            coords,
                                            direction = "RefToFocal",
                                            AAinput = FALSE,
                                            AAoutput = FALSE)
    Nt_mutation <- mutations_list$Nt_mutation[j]
    if (!is.na(Nt_pos)) {  # only proceed if corresponding Nt presents in target sequence!
      # result$AA_target[i] <- as.character(protein_target[AA_pos])
      result$Nt_target[j] <- as.character(target_sequence[Nt_pos])
      # nt_pos <- AA_pos * 3 - 2
      # result$Codon_target[i] <- as.character(target_sequence[nt_pos:(nt_pos + 2)])
      # result$n_possible[i] <- n_mutations_in_codon(AA_pos, 
      #                                              mutations_list$AA_mutation[i], 
      #                                              target_sequence)
      result$n_possible[j] <- 1
    }
  }
  result$alig_score <- pwalign::score(coords_alig$alignment) # alignment score for each target sequence
  
  # calculate Levenshtein distance between E. coli and target gene core region:
  # determine whether to use the core region or the original sequence
  
  from_target <- ALJEbinf::translateCoordinate(510, coords, direction = "RefToFocal")
  to_target <- ALJEbinf::translateCoordinate(920, coords, direction = "RefToFocal")
  core_target <- target_sequence[from_target:to_target]
  core_Ecoli <- reference_Ecoli[510:920]
  result$core_dist <- pwalign::stringDist(c(Biostrings::DNAStringSet(core_target), Biostrings::DNAStringSet(core_Ecoli)),
                                          method = "levenshtein")
  return(result)
}

#' Function to process raw output from the mutation screen
#'
#' @param output raw mutation screen output.
#'
#' @return A table with optimised species names and new columns "genus",
#' "mutation_name" and "mutation_category".
#' @export
#'
process_output_nt <- function(output) {
  processed_output <- output |>
    mutate(species= gsub("\\[|\\]", "", species)) |> # remove [] sign from species names
    mutate(species = gsub("\\'|\\'", "", species)) |> # remove '' sign from species names
    mutate(species = str_replace_all(species, "Candidatus_", "")) |> # remove this name from all species names
    mutate(species = gsub("_", " ", species)) |> # remove underscore from species names
    mutate(mutation_name = paste(Nt_pos_Ecoli, Nt_mutation, sep = "_")) |> # create a unique name for mutations
    mutate(genus = sub(" .*", "", species)) |>
    #categorise each mutation in three possible options:
    # 1) mutation is already present ("present")
    # 2) mutation is not present but can be gained ("possible")
    # 3) mutation is not present and also cannot be gained ("impossible")
    mutate(mutation_category = ifelse(Nt_mutation == Nt_target, 
                                      "present",
                                      ifelse(n_possible == 0L,
                                             "impossible",
                                             "possible")))
  return(processed_output)
}


#' filters final_output_nt to retrieve data related to high quality sequences and reliable mutation list
#' @param output a data frame providing the results of mutation screening in all extracted target sequences
#' @param min_seq_length a numeric vector providing the minimum length of an acceptable target sequence
#' @param max_dist_from_ref a numeric vector providing the maximum distance of an acceptable target sequence against reference sequence
#'
#' @return a data frame providing the results of mutation screening of high quality target sequences across bacterial species
#' @export
#'
#' @examples filter_output_nt(output, 2000, -Inf, 1500)
filter_output_nt <- function(output, min_seq_length, min_alig_score, max_core_dist) {
  filtered_output <- output |>
    filter(target_length >= min_seq_length) |>
    filter(alig_score >= min_alig_score) |> 
    filter(core_dist <= max_core_dist) |> 
    filter(!is.na(Nt_target), !is.na(n_possible)) |> 
    unique() 
  return(filtered_output)
}
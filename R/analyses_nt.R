#' Prepare a (possibly filtered) table of mutations to be screened:
#'
#' @param muts a data frame providing info on rrs (nucleotide) mutations reported in the literature
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
    #categorise each mutation in two possible options:
    # 1) mutation is already present ("present")
    # 2) mutation is not present ("possible")
    # mutate(mutation_category = ifelse(Nt_mutation == Nt_target, 
    #                                   "present",
    #                                   ifelse(n_possible == 0L,
    #                                          "impossible",
    #                                          "possible")))
    mutate(mutation_category = ifelse(Nt_mutation == Nt_target, 
                                      "present", "possible"))
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


#' Compare multiple rrs gene copies within bacterial species
#'
#' This function identifies bacterial species that carry more than one copy of the rrs gene 
#' and performs pairwise comparisons of all copies to assess sequence similarity and core-region divergence.
#' It calculates alignment scores and Levenshtein distances for the E. coli 16S rRNA core region (nt 510–920)
#' across all copy pairs, and summarizes resistance status and pairwise metrics for each species.
#'
#' @param filtered_output A data frame containing high-quality rrs sequence metadata, including accession numbers, species, genus, gene copy IDs, target names, and mutation information.
#' @param rrs_target_sequences A named list of DNA sequences corresponding to each rrs target in filtered_output.
#' @param rrs_reference_Ecoli A reference E. coli rrs sequence used to define the core region coordinates for pairwise comparisons.
#' @param outfile Character string specifying the path to the CSV file where the summarized results will be written.
#' @export
#' 

compare_rrs_copies <- function(filtered_output, rrs_target_sequences, rrs_reference_Ecoli) {
  
  # Create file with header if it doesn't exist
  if (!file.exists(outfile)) {
    write_csv(tibble(
      accession_number = character(),
      species_name = character(),
      genus = character(),
      resistance_status = integer(),
      mean_alig_score = double(),
      mean_core_dist = double(),
      n_pairs = integer()
    ), outfile)
  }
  
  # Identify species with more than one rrs copy
  multicopy_species <- filtered_output |>
    distinct(accession_numbers, gene_copy) |>
    count(accession_numbers, name = "n") |>
    filter(n > 1L) |>
    pull(accession_numbers)
   
  # Define E. coli core region coordinates (510–920 nt)
  core_Ecoli <- rrs_reference_Ecoli[510:920]
  
  plan(multisession, workers = 5)  # parallel execution
  
  multiseq_stats <- future_lapply(seq_along(multicopy_species), function(i) {
    acc <- multicopy_species[i]
    cat(sprintf("Working on multi-copy species #%d/%d...\n", i, length(multicopy_species)))
    # Extract metadata
    seq_names <- filtered_output |>
      filter(accession_numbers == acc) |>
      pull(target_name) |>
      unique()
    
    species_name <- filtered_output |>
      filter(accession_numbers == acc) |>
      pull(species) |>
      unique()
    
    genus <- filtered_output |>
      filter(accession_numbers == acc) |>
      pull(genus) |>
      unique()
    
    # Determine resistance status (how many copies carry resistance mutations)
    resistance_status <- filtered_output |>
      filter(accession_numbers == acc) |>
      group_by(gene_copy) |>
      summarise(resistant_copy = any(mutation_category == "present"), .groups = "drop") |>
      pull(resistant_copy) |>
      sum()
    
    # All pairwise combinations of copies
    combs <- combn(seq_names, 2, simplify = FALSE)
    
    results <- lapply(combs, function(pair) {
      seq1 <- rrs_target_sequences[[pair[1]]]
      seq2 <- rrs_target_sequences[[pair[2]]]
      
      # Map E. coli core region onto the two sequences
      coords1 <- getCoordinatesNt(seq1, rrs_reference_Ecoli)
      coords2 <- getCoordinatesNt(seq2, rrs_reference_Ecoli)
      
      from1 <- ALJEbinf::translateCoordinate(510, coords1, direction = "RefToFocal")
      to1   <- ALJEbinf::translateCoordinate(920, coords1, direction = "RefToFocal")
      core1 <- seq1[from1:to1]
      
      from2 <- ALJEbinf::translateCoordinate(510, coords2, direction = "RefToFocal")
      to2   <- ALJEbinf::translateCoordinate(920, coords2, direction = "RefToFocal")
      core2 <- seq2[from2:to2]
      
      tibble(
        accession_number = acc,
        alig_score = pwalign::score(
          pwalign::pairwiseAlignment(seq1, seq2)
        ),
        core_dist = as.double(
          pwalign::stringDist(
            c(Biostrings::DNAStringSet(core1), Biostrings::DNAStringSet(core2)),
            method = "levenshtein"
          )
        )
      )
    }) |> bind_rows()
    
    # Summarise pairwise results for this species
    summary <- results |>
      summarise(
        accession_number = first(accession_number),
        species_name     = first(species_name),
        genus            = first(genus),
        resistance_status = resistance_status,
        mean_alig_score  = mean(alig_score, na.rm = TRUE),
        mean_core_dist   = mean(core_dist, na.rm = TRUE),
        n_pairs          = n()
      )
    
    # Append result to file immediately
    write_csv(summary, outfile, append = TRUE)
    return(summary)
  })
  
  # Combine all results
  output <- bind_rows(multiseq_stats)
  return(multiseq_stats)
}

#' Summarise output to species-level
#'
#' @param output 
#'
#' @return A data frame (tibble) with one row per species, containing columns for 
#' instrinsic resistance ("resistance"), evolvabilityI ("evolvabilityI", defined as the
#' number of mutations that are accessible).
#' @export
#'
get_species_output_nt <- function(output) {
  species_output <-  output |>
    # step 1: merge multiple gene copies into one:
    group_by(species, genus, accession_numbers, mutation_name) |>
    summarise(mutation_category = ifelse(any(mutation_category == "present"),
                                         "present",
                                         ifelse(any(mutation_category == "possible"), "possible", "impossible")),
              n_possible = sum(n_possible),
              .groups = "drop") |>
    # step 2: summarise across all screened mutations:
    group_by(species, genus, accession_numbers) |>
    summarise(resistance = any(mutation_category == "present"),
              # evolvabilityI = sum(mutation_category == "possible"),
              # evolvabilityII = sum(n_possible),
              .groups = "drop") |>
    # step 3: discard sequences in cases where there are >1 genome per species:
    group_by(species, genus) |>
    summarise(resistance = dplyr::first(resistance),
              # evolvabilityI = dplyr::first(evolvabilityI),
              # evolvabilityII = dplyr::first(evolvabilityII),
              .groups = "drop") |>
    # step 4: turn resistance column into categorical:
    mutate(resistance = ifelse(resistance, "resistant", "susceptible"))
  return(species_output)
}

#' Calculate conservation and distance scores along the E. coli gene sequence
#' 
#' For each position within the gene sequence, four scores are calculated for each gene sequence:
#' 
#' hamming_Ecoli: whether or not the nucleotide is different from the one in E. coli
#' hamming_rnd: whether or not the nucleotide is different between two randomly chosen gene sequences
#'
#' @param target_sequences 
#' @param reference_Ecoli 
#' @param alig_path 
#'
#' @returns A list with four matrices (one for each of the scores described above) and 
#' a tibble containing all the means for each position, across sequences
#' 
get_conservation_nt <- function(target_sequences, 
                             reference_Ecoli,
                             n_rnd = 1e5,
                             alig_path = "./output/alignments_nt",
                             n_workers)
{
  # l <- length(translate(reference_Ecoli))
  l <- length(reference_Ecoli)
  m <- length(target_sequences)
  # reference_Ecoli_AA <- translate(reference_Ecoli)
  # grantham <- granthamMatrix()
  
  plan(multisession, workers = n_workers) # to parallelise the function
  
  # distances to E. coli:
  results_Ecoli <- future_lapply(1:m, function(k) {
    hamming_Ecoli <- rep(NA, l)
    # grantham_Ecoli <- rep(NA, l)
    try( {
      load(paste0(alig_path, "/", names(target_sequences)[k], ".RData"))
      coords <- coords_alig$coordinates
      rrs_target_nt <- target_sequences[[k]]
      
      for(pos_Ecoli in 1:(l-1)) { # loop through E. coli rpsL
        pos_target <- translateCoordinate(pos_Ecoli, coords, direction = "RefToFocal")
        if (!is.na(pos_target)) {
          nt_target <- as.character(rrs_target_nt[pos_target])
          nt_Ecoli <- as.character(reference_Ecoli[pos_Ecoli])
          hamming_Ecoli[pos_Ecoli] <- (nt_target != nt_Ecoli)
          # if ((nt_Ecoli %in% colnames(grantham)) && (nt_target %in% colnames(grantham))) {
          #   grantham_Ecoli[pos_Ecoli] <- grantham[nt_target, nt_Ecoli]
          # } else {
          #   warning(paste0("Unknown nucleotide (", nt_target, ")."))
          # }
        } else {
          hamming_Ecoli[pos_Ecoli] <- TRUE
          # grantham_Ecoli[pos_Ecoli] <- NA
        }
      }
    }
    )
    return(list(hamming_Ecoli = hamming_Ecoli))
  }
  )
  
  # distances between random sequences:
  ijs <- expand_grid(i = 1:m, j = 1:m) |>
    filter(i < j) |>
    slice_sample(n = n_rnd)
  
  results_rnd <- future_lapply(1:n_rnd, function(k) {
    hamming_rnd <- rep(NA, l)
    # grantham_rnd <- rep(NA, l)
    try( {
      # sequence 1:
      i <- ijs$i[k]
      load(paste0(alig_path, "/", names(target_sequences)[i], ".RData"))
      coords1 <- coords_alig$coordinates
      rrs_target1_nt <- target_sequences[[i]]
      
      # sequence 2:
      j <- ijs$j[k]
      load(paste0(alig_path, "/", names(target_sequences)[j], ".RData"))
      coords2 <- coords_alig$coordinates
      rrs_target2_nt <- target_sequences[[j]]
      
      for(pos_Ecoli in 1:(l-1)) { # loop through E. coli rpsL
        pos_target1 <- translateCoordinate(pos_Ecoli, coords1, direction = "RefToFocal")
        pos_target2 <- translateCoordinate(pos_Ecoli, coords2, direction = "RefToFocal")
        if ((!is.na(pos_target1)) && (!is.na(pos_target2))) {
          nt_target1 <- as.character(rrs_target1_nt[pos_target1])
          nt_target2 <- as.character(rrs_target2_nt[pos_target2])
          hamming_rnd[pos_Ecoli] <- (nt_target1 != nt_target2)
          # if ((AA_target1 %in% colnames(grantham)) && (AA_target2 %in% colnames(grantham))) {
          #   grantham_rnd[pos_Ecoli] <- grantham[AA_target1, AA_target2]
          # } else {
          #   warning(paste0("Unknown AA (", AA_target, ")."))
          # }
        } else if (xor(is.na(pos_target1), is.na(pos_target2))) {
          hamming_rnd[pos_Ecoli] <- TRUE
          # grantham_rnd[pos_Ecoli] <- NA
        } else {
          hamming_rnd[pos_Ecoli] <- NA
          # grantham_rnd[pos_Ecoli] <- NA
        }
      }
    }
    )
    return(list(hamming_rnd = hamming_rnd))
  }
  )
  
  hamming_Ecoli <- matrix(NA, nrow = m, ncol = l,
                          dimnames = list(names(target_sequences), paste0("P", 1:l)))
  # grantham_Ecoli <- hamming_Ecoli
  for(i in 1:m) {
    hamming_Ecoli[i,] <- results_Ecoli[[i]]$hamming_Ecoli
    # grantham_Ecoli[i,] <- results_Ecoli[[i]]$grantham_Ecoli
  }
  
  hamming_rnd <- matrix(NA, nrow = n_rnd, ncol = l,
                        dimnames = list(paste0(names(target_sequences)[ijs$i], 
                                               ":::",
                                               names(target_sequences)[ijs$j]),
                                        paste0("P", 1:l)))
  # grantham_rnd <- hamming_rnd
  for(i in 1:n_rnd) {
    hamming_rnd[i,] <- results_rnd[[i]]$hamming_rnd
    # grantham_rnd[i,] <- results_rnd[[i]]$grantham_rnd
  }
  
  # means aross all sequences:
  means = tibble(hamming_Ecoli = colMeans(hamming_Ecoli, na.rm = TRUE),
                 # grantham_Ecoli = colMeans(grantham_Ecoli, na.rm = TRUE),
                 hamming_rnd = colMeans(hamming_rnd, na.rm = TRUE))
                 # grantham_rnd = colMeans(grantham_rnd, na.rm = TRUE))
  
  return(list(hamming_Ecoli = hamming_Ecoli, 
              # grantham_Ecoli = grantham_Ecoli,
              hamming_rnd = hamming_rnd, 
              # grantham_rnd = grantham_rnd,
              means = means))
}






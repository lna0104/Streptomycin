#' Number of mutations translated to same amino acid in a codon
#'
#' @param AA_pos   a value indicating the amino acid position within a DNA sequence
#' @param mut_AA   a substituted amino acid
#' @param DNA_seq  a DNA string object
#'
#' @return the value indicating all possible mutations which bring the same amino acid substitution within a protein sequence
#' @export
#'
#' @examples n_mutations_in_codon(2, "D", DNA_string)
n_mutations_in_codon <- function(AA_pos, mut_AA, DNA_seq) {
  codon <- substr(DNA_seq, AA_pos*3 - 2, AA_pos*3)
  mut_codons <- mutant_codons(as.character(codon))
  n_muts <- sum(mut_codons$AA == mut_AA)
  return(n_muts)
}


#' Possible amino acids result from all possible point mutations in any three-letter codon 
#'
#' @param codon a three-letter DNA string object which can be translated to an amino acid
#'
#' @return a data frame of all nine possible codons and amino acids for any three-letter codon  
#' @export
#'
#' @examples mutant_codon(TCG)
mutant_codons <- function(codon){
  
  ## check that argument codon is a three-letter string of A, C, T, G:
  
  if (str_length(codon) != 3){
    return("string is of wrong length")
  }
  
  ## make empty data frame
  mutants <- data.frame(codon = rep(NA, 9),
                        AA = NA)
  
  ## loop through all positions to get all nine possible nucleotide composition in set of three nuleotides
  nts <- c("A", "T", "G", "C")
  m <- 1
  for(p in 1:3) 
    for(n in 1:4) {
      if (str_sub(codon, p, p) != nts[n]) {
        mutant_codon <- codon
        # codon
        str_sub(mutant_codon, p, p) <- nts[n]
        mutants$codon[m] <- mutant_codon
        # add AA here:
        mutants$AA[m] <- as.character(translate(DNAString(mutant_codon), no.init.codon = TRUE))
        m <- m + 1
        
      }
    }
  return(mutants)
}



#' Title: check out mutation
#'
#' @param mutations_list a data frame providing all amino acid substitutions reported in Google sheet based on coordinations
#' @param target_sequences a list of multiple DNAString objects (gene sequences retrieved from different bacterial species) 
#' @param reference_Ecoli a DNAString object (gene sequence of Escherichia coli MG1655)
#' @param target_gene a name of target gene
#' @return A data frame of all amino acid substitutions within all target sequences and check if there is a amino acid substantiation and whether amino acid substantiation is new or has been reported before or there is a different amino acid in the original position 
#' @export
#'
#' @examples check_out_mutation(mutations_list, target_sequence, reference_Ecoli, target_gene)
check_out_mutation <- function(mutations_list, 
                               target_sequence, 
                               reference_Ecoli,
                               target_gene,
                               alig_file = NULL) {
  if (is.null(target_sequence) || length(target_sequence) == 0) {
    stop(sprintf('Empty "%s" sequence.', target_gene))
  }
  coords_alig <- ALJEbinf::getCoordinates(target_sequence, reference_Ecoli, aligOutput = TRUE)
  if (!is.null(alig_file)) {
    save(coords_alig, file = alig_file)
  }
  coords <- coords_alig$coordinates
  protein_target <- Biostrings::translate(target_sequence)
  protein_Ecoli <- Biostrings::translate(reference_Ecoli)
  
  result <- mutations_list
  result$AA_target <- NA # what is the amino acid at that position in focal species
  result$Codon_target <- NA # what's the codon at that position in focal species
  result$n_possible <- NA # number of single nt mutations in sequence that would produce the AA mutation
  for (i in 1:nrow(mutations_list)) {
    AA_pos <- ALJEbinf::translateCoordinate(mutations_list$AA_pos_Ecoli[i], 
                                            coords,
                                            direction = "RefToFocal",
                                            AAinput = TRUE,
                                            AAoutput = TRUE)
    AA_mutation <- mutations_list$AA_mutation[i]
    if (!is.na(AA_pos)) {  # only proceed if corresponding AA presents in target sequence!
      result$AA_target[i] <- as.character(protein_target[AA_pos])
      nt_pos <- AA_pos * 3 - 2
      result$Codon_target[i] <- as.character(target_sequence[nt_pos:(nt_pos + 2)])
      result$n_possible[i] <- n_mutations_in_codon(AA_pos, 
                                                   mutations_list$AA_mutation[i], 
                                                   target_sequence)
    }
  }
  result$alig_score <- pwalign::score(coords_alig$alignment) # alignment score for each target sequence
  
  # calculate Levenshtein distance between E. coli and target gene core region:
  # determine whether to use the core region or the original sequence
  if (length(protein_Ecoli) > 300) {
    # Use core region if gene length > 300
    from_target <- ALJEbinf::translateCoordinate(501, coords, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
    to_target <- ALJEbinf::translateCoordinate(600, coords, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
    core_target <- protein_target[from_target:to_target]
    core_Ecoli <- protein_Ecoli[501:600]
  } else {
    # Use full original sequence if gene length ≤ 300
    core_target <- protein_target
    core_Ecoli <- protein_Ecoli
  }
  result$core_dist <- pwalign::stringDist(c(Biostrings::AAStringSet(core_target),Biostrings::AAStringSet(core_Ecoli)),
                                             method = "levenshtein")
  return(result)
}

#' Screen a set of target sequences for existence and evolvability over list of mutations
#'
#' @param target_sequences A DNAStringSet object containing all target gene sequences
#' @param mutation_list A data frame containing all the mutations to be screened for
#'
#' @return A list of either data frames (if screen for a sequence was successful) or error messages (if not).
#'
#'
screen_target_sequences <- function(target_sequences, reference_Ecoli, mutation_list, target_gene, n_workers) {
  # Ensure required directories exist
  if (!dir.exists("./output/alignments")) dir.create("./output/alignments", recursive = TRUE)
  if (!dir.exists("./output/progress")) dir.create("./output/progress", recursive = TRUE)
  
  plan(multisession, workers = n_workers) # to parallelise the function
  final_output <- future_lapply(1:length(target_sequences), function(i) {
    cat(paste0("Working on species #", i, "/", length(target_sequences), "...\n")) #notify when the analysis of each sequence starts
    output <- tryCatch(
      check_out_mutation(mutation_list,              # compare target sequences with reference sequence to identify the amino acid at each position
                         target_sequences[[i]], 
                         reference_Ecoli,
                         target_gene,
                         alig_file = paste0("./output/alignments/", names(target_sequences)[[i]], ".RData")),
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
    
    saveRDS(i, file = paste0("./output/progress/sequence_", i, "_out_of_", length(rpsL_target_sequences), "_complete.rds"))
    return(output)
  })
  return(final_output)
}


#' Prepare a (possibly filtered) table of mutations to be screened:
#'
#' @param muts a data frame providing info on rpoB mutations reported in the literature
#' @param min_n_species a numeric vector providing the minimum number of species which a mutation has been reported in
#' @param min_n_studies a numeric vector providing the minimum number studies which a mutation has been reported in
#' @param origin Origin of mutations to be included, e.g., "Isolate" or "Lab mutant". 
#' The default value of "ANY" means that all origins will be included. 
#'
#' @return a data frame providing info on filtered reported mutations in the literature(reliable mutations)
#' @export
#'
#' @examples filter_mutations(muts, min_n_species = 5)
filter_mutations <- function(muts, 
                             min_n_species = 1, 
                             min_n_studies = 1, 
                             origin = "ANY") {
  
  mutation_list <- muts |>
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_original), !is.na(AA_mutation)) |>
    filter(is.na(Warning) | Warning == "AA_pos inconsistent with AA_pos_Ecoli")
  
  if (origin == "ANY") {
    origin <- muts |> pull(Origin) |> unique()
  } 
  
  mutation_list_species <- mutation_list |>
    select(Species, AA_pos_Ecoli, AA_mutation) |>
    distinct() |>
    group_by(AA_pos_Ecoli, AA_mutation) |>
    summarise(n_species = n(), .groups = "drop")
  
  mutation_list_studies <- mutation_list |>
    select(Ref_code, AA_pos_Ecoli, AA_mutation) |>
    distinct() |>
    group_by(AA_pos_Ecoli, AA_mutation) |>
    summarise(n_studies = n(), .groups = "drop")
  
  mutation_list <- mutation_list_studies |>
    full_join(mutation_list_species, by = join_by(AA_pos_Ecoli, AA_mutation)) |>
    filter(n_studies >= min_n_studies) |>
    filter(n_species >= min_n_species) |>
    select(AA_pos_Ecoli, AA_mutation)
  
  return(mutation_list)
}

#' Function to process raw output from the mutation screen
#'
#' @param output raw mutation screen output.
#'
#' @return A table with optimised species names and new columns "genus",
#' "mutation_name" and "mutation_category".
#' @export
#'
process_output <- function(output) {
  processed_output <- output |>
    mutate(species= gsub("\\[|\\]", "", species)) |> # remove [] sign from species names
    mutate(species = gsub("\\'|\\'", "", species)) |> # remove '' sign from species names
    mutate(species = str_replace_all(species, "Candidatus_", "")) |> # remove this name from all species names
    mutate(species = gsub("_", " ", species)) |> # remove underscore from species names
    mutate(mutation_name = paste(AA_pos_Ecoli, AA_mutation, sep = "_")) |> # create a unique name for mutations
    mutate(genus = sub(" .*", "", species)) |>
    #categorise each mutation in three possible options:
    # 1) mutation is already present ("present")
    # 2) mutation is not present but can be gained ("possible")
    # 3) mutation is not present and also cannot be gained ("impossible")
    mutate(mutation_category = ifelse(AA_mutation == AA_target, 
                                      "present",
                                      ifelse(n_possible == 0L,
                                             "impossible",
                                             "possible")))
  return(processed_output)
}


#' filters final_output to retrieve data related to high quality sequences and reliable mutation list
#' @param output a data frame providing the results of mutation screening in all extracted target sequences
#' @param min_seq_length a numeric vector providing the minimum length of an acceptable target sequence
#' @param max_dist_from_ref a numeric vector providing the maximum distance of an acceptable target sequence against reference sequence
#'
#' @return a data frame providing the results of mutation screening of high quality target sequences across bacterial species
#' @export
#'
#' @examples filter_output(output, 2000, -Inf, 1500)
filter_output <- function(output, min_seq_length, min_alig_score, max_core_dist) {
  filtered_output <- output |>
    filter(target_length >= min_seq_length) |>
    filter(alig_score >= min_alig_score) |> 
    filter(core_dist <= max_core_dist) |> 
    filter(!is.na(AA_target), !is.na(Codon_target), !is.na(n_possible)) |> 
    unique() 
  return(filtered_output)
}



compare_gene_copies <- function(filtered_output, target_sequences, reference_Ecoli) {

  # species with more than one gene copy:
  multicopy_species <- filtered_output |>
    select(accession_numbers, gene_copy) |>
    distinct() |>
    group_by(accession_numbers) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n > 1L) |>
    pull(accession_numbers)
  
  # If there are no species with multiple copies, proceed to the next step
  if (length(multicopy_species) == 0) {
    message("No species with multiple gene copies found.")
    return(invisible(NULL))  
  }
  
  protein_Ecoli <- Biostrings::translate(reference_Ecoli)
  # If the gene is rpoB
  if (length(protein_Ecoli) > 300) {
    core_Ecoli <- protein_Ecoli[501:600]
  }else{
  # If the gene is rpsL
    core_Ecoli <- protein_Ecoli
  }
  
  plan(multisession) # to parallelise the function
  multiseq_stats <- future_lapply(1:length(multicopy_species), function(i) {
    cat(paste0("Working on multi-copy species #", i, "/", length(multicopy_species), "...\n")) #notify when the analysis of each sequence starts
    seq_names <- filtered_output |>
      filter(accession_numbers == multicopy_species[i]) |>
      pull(target_name) |>
      unique()
    species_name <- filtered_output |>
      filter(accession_numbers == multicopy_species[i]) |>
      pull(species) |>
      unique()
    genus <- filtered_output |>
      filter(accession_numbers == multicopy_species[i]) |>
      pull(genus) |>
      unique()
    # determine resistance status (how many of the two gene copies confer resistance)
    resistance_status <- filtered_output |>
      filter(accession_numbers == multicopy_species[i]) |>
      group_by(gene_copy) |>
      summarise(resistant_copy = any(mutation_category == "present"), .groups = "drop") |>
      pull(resistant_copy) |>
      sum()  
    
    if (length(seq_names) > 2L) {
      warning("More than two gene copies in genome ", multicopy_species[i], ", ignored.")
      output <- NULL
    } else {
      coords1 <- ALJEbinf::getCoordinates(target_sequences[[seq_names[1]]], 
                                          reference_Ecoli)
      coords2 <- ALJEbinf::getCoordinates(target_sequences[[seq_names[2]]], 
                                          reference_Ecoli)
      protein_target1 <- Biostrings::translate(target_sequences[[seq_names[1]]])
      protein_target2 <- Biostrings::translate(target_sequences[[seq_names[2]]])
      
      if (length(protein_target1) > 300 & length(protein_target2) > 300) {
        # calculate Levenshtein distance between E. coli and target rpoB core region (AA pos 501...600):
        from_target1 <- ALJEbinf::translateCoordinate(501, coords1, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
        to_target1 <- ALJEbinf::translateCoordinate(600, coords1, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
        core_target1 <- protein_target1[from_target1:to_target1]
        from_target2 <- ALJEbinf::translateCoordinate(501, coords2, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
        to_target2 <- ALJEbinf::translateCoordinate(600, coords2, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
        core_target2 <- protein_target2[from_target2:to_target2]
      }else{
        # calculate Levenshtein distance between E. coli and target rpsL
        core_target1 <- protein_target1
        core_target2 <- protein_target2
      }
      
      
      output <- list(accession_number = multicopy_species[i],
                     species_name = species_name,
                     genus = genus,
                     resistance_status = resistance_status,
                     alig_score = pwalign::score(pwalign::pairwiseAlignment(protein_target1, protein_target2)),
                     core_dist = as.double(pwalign::stringDist(c(Biostrings::AAStringSet(core_target1), 
                                                                 Biostrings::AAStringSet(core_target2)),
                                                               method = "levenshtein")))
    }
    return(output)
  })
  return(do.call(rbind, purrr::map(multiseq_stats, as_tibble)))
}

#' Summarise output to species-level
#'
#' @param output 
#'
#' @return A data frame (tibble) with one row per species, containing columns for 
#' instrinsic resistance ("resistance"), evolvability I ("evolvabilityI", defined as the
#' number of AA mutations that are accessible), and evolvability II ("evolvabilityII", 
#' defined as the total number of mutations that can produce a resistance AA mutation).
#' @export
#'
get_species_output <- function(output) {
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
              evolvabilityI = sum(mutation_category == "possible"),
              evolvabilityII = sum(n_possible),
              .groups = "drop") |>
    # step 3: discard sequences in cases where there are >1 genome per species:
    group_by(species, genus) |>
    summarise(resistance = dplyr::first(resistance),
              evolvabilityI = dplyr::first(evolvabilityI),
              evolvabilityII = dplyr::first(evolvabilityII),
              .groups = "drop") |>
    # step 4: turn resistance column into categorical:
    mutate(resistance = ifelse(resistance, "resistant", "susceptible"))
  return(species_output)
}

# The following function returns the minimum and maximum possible
# values of both evolvability I and II
get_theoretical_evolvabilities <- function(output) {
  mutation_list <- output |>
    select(AA_pos_Ecoli, AA_mutation) |>
    distinct() |>
    mutate(mut_name = paste(AA_pos_Ecoli, AA_mutation, sep = '_'))
  AA_pos_Ecoli <- mutation_list |>
    pull(AA_pos_Ecoli) |>
    unique()
  nts <- c('A', 'C', 'T', 'G')
  all_codons <- paste0(rep(nts, each = 16),
                       rep(paste0(rep(nts, each = 4), rep(nts, 4)), 4))
  possibilities <- expand_grid(AA_pos_Ecoli, 
                                codon_from = all_codons, 
                                codon_to = all_codons) |>
    mutate(hamming = hamming_str(codon_from, codon_to)) |>
    filter(hamming == 1L) |>
    select(-hamming) |>
    mutate(AA_from = Biostrings::GENETIC_CODE[codon_from],
           AA_to = Biostrings::GENETIC_CODE[codon_to]) |>
    mutate(AA_mut_name_from = paste(AA_pos_Ecoli, AA_from, sep = '_'),
           AA_mut_name_to = paste(AA_pos_Ecoli, AA_to, sep = '_')) |>
    mutate(resistant_from = AA_mut_name_from %in% mutation_list$mut_name,
           resistant_to = AA_mut_name_to %in% mutation_list$mut_name) |>
    filter(AA_from != '*')
  
  evolvabilitiesI <- possibilities |>
    filter(!resistant_from) |>
    select(AA_pos_Ecoli, codon_from, AA_from, AA_mut_name_from, AA_mut_name_to, resistant_to) |>
    distinct() |>
    group_by(AA_pos_Ecoli, codon_from, AA_from, AA_mut_name_from) |>
    summarise(evolvabilityI = sum(resistant_to), .groups = 'drop') |>
    group_by(AA_pos_Ecoli) |>
    summarise(min_evolvabilityI = min(evolvabilityI),
              max_evolvabilityI = max(evolvabilityI))
  
  evolvabilitiesI <- possibilities |>
    filter(!resistant_from) |>
    select(AA_pos_Ecoli, codon_from, AA_from, AA_mut_name_from, AA_mut_name_to, resistant_to) |>
    distinct() |>
    group_by(AA_pos_Ecoli, codon_from, AA_from, AA_mut_name_from) |>
    summarise(evolvabilityI = sum(resistant_to), .groups = 'drop') |>
    group_by(AA_pos_Ecoli) |>
    summarise(min_evolvabilityI = min(evolvabilityI),
              max_evolvabilityI = max(evolvabilityI))
  evolvabilitiesII <- possibilities |>
    group_by(AA_pos_Ecoli, codon_from, AA_from, AA_mut_name_from) |>
    summarise(evolvabilityII = sum(resistant_to), .groups = 'drop') |>
    group_by(AA_pos_Ecoli) |>
    summarise(min_evolvabilityII = min(evolvabilityII),
              max_evolvabilityII = max(evolvabilityII))
  
  return(list(
    min_evolvabilityI = evolvabilitiesI |> pull(min_evolvabilityI) |> sum(),
    max_evolvabilityI = evolvabilitiesI |> pull(max_evolvabilityI) |> sum(),
    min_evolvabilityII = evolvabilitiesII |> pull(min_evolvabilityII) |> sum(),
    max_evolvabilityII = evolvabilitiesII |> pull(max_evolvabilityII) |> sum()
    ))
}


get_resistance_taxonomy <- function(output, bacterial_taxonomy, file_path) {
  resistant_species <- output |>
    group_by(genus, species, accession_numbers) |>
    summarise(resistant = any(mutation_category == "present"), .groups = "drop") |>
    left_join(bacterial_taxonomy, join_by(genus == genus)) 
  
  resistant_genera <- resistant_species |>
    group_by(genus, family, order, class, phylum) |>
    summarise(n = n(), n_res = sum(resistant), .groups = "drop") |>
    mutate(f_res = n_res/n)
  
  resistant_families <- resistant_species |>
    group_by(family, order, class, phylum) |>
    summarise(n = n(), n_res = sum(resistant), .groups = "drop") |>
    mutate(f_res = n_res/n)
  
  resistant_orders <- resistant_species |>
    group_by(order, class, phylum) |>
    summarise(n = n(), n_res = sum(resistant), .groups = "drop") |>
    mutate(f_res = n_res/n)
  
  resistant_classes <- resistant_species |>
    group_by(class, phylum) |>
    summarise(n = n(), n_res = sum(resistant), .groups = "drop") |>
    mutate(f_res = n_res/n)
  
  write_csv(resistant_genera, paste0(file_path, "resistant_genera.csv"))
  write_csv(resistant_families, paste0(file_path, "resistant_families.csv"))
  write_csv(resistant_orders, paste0(file_path, "resistant_orders.csv"))
  write_csv(resistant_classes, paste0(file_path, "resistant_classes.csv"))
}


#' Find mismatches between wildtype AAs in reports vs. aligned sequences from NCBI
#'
#' @param muts Completed table of reported mutations
#' @param output Ouput from the mutation screen
#'
#' @return A table with all mismatches
#' @export
#'
get_wt_AA_mismatches <- function(muts, output) {
  output_reduced <- output |>
    select(species, AA_pos_Ecoli, AA_mutation, AA_target, Codon_target)
  mismatches <- muts |>
    left_join(output_reduced, 
              by = join_by(Species == species, AA_pos_Ecoli == AA_pos_Ecoli, AA_mutation == AA_mutation), 
              relationship = "many-to-many") |>
    select(Species, AA_pos_Ecoli, AA_mutation, AA_original, AA_target, Ref_code, Comments) |>
    filter(AA_original != AA_target)
  return(mismatches)
}


#' Calculate conservation and distance scores along the E. coli gene sequence
#' 
#' For each position within the gene sequence, four scores are calculated for each gene sequence:
#' 
#' hamming_Ecoli: whether or not the AA is different from the one in E. coli
#' hamming_rnd: whether or not the AA is different between two randomly chosen gene sequences
#' grantham_Ecoli: distance between AA and AA in E. coli, according to the Grantham distance matrix
#' grantham_rnd: distance between AAs in two randomly chosen gene sequences, according to the Grantham distance matrix
#'
#' @param target_sequences 
#' @param reference_Ecoli 
#' @param alig_path 
#'
#' @returns A list with four matrices (one for each of the scores described above) and 
#' a tibble containing all the means for each position, across sequences
#' 
get_conservation <- function(target_sequences, 
                             reference_Ecoli,
                             n_rnd = 1e5,
                             alig_path = "./output/alignments")
{
  l <- length(translate(reference_Ecoli))
  m <- length(target_sequences)
  reference_Ecoli_AA <- translate(reference_Ecoli)
  grantham <- granthamMatrix()
  
  plan(multisession, workers = 10) # to parallelise the function
  
  # distances to E. coli:
  results_Ecoli <- future_lapply(1:m, function(k) {
    hamming_Ecoli <- rep(NA, l)
    grantham_Ecoli <- rep(NA, l)
    try( {
      load(paste0(alig_path, "/", names(target_sequences)[k], ".RData"))
      coords <- coords_alig$coordinates
      rpoB_target_AA <- translate(target_sequences[[k]])

      for(pos_Ecoli in 1:(l-1)) { # loop through E. coli rpoB
        pos_target <- translateCoordinate(pos_Ecoli, coords, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
        if (!is.na(pos_target)) {
          AA_target <- as.character(rpoB_target_AA[pos_target])
          AA_Ecoli <- as.character(reference_Ecoli_AA[pos_Ecoli])
          hamming_Ecoli[pos_Ecoli] <- (AA_target != AA_Ecoli)
          if ((AA_Ecoli %in% colnames(grantham)) && (AA_target %in% colnames(grantham))) {
            grantham_Ecoli[pos_Ecoli] <- grantham[AA_target, AA_Ecoli]
          } else {
            warning(paste0("Unknown AA (", AA_target, ")."))
          }
        } else {
          hamming_Ecoli[pos_Ecoli] <- TRUE
          grantham_Ecoli[pos_Ecoli] <- NA
        }
      }
    }
    )
    return(list(hamming_Ecoli = hamming_Ecoli, grantham_Ecoli = grantham_Ecoli))
  }
  )
  
  # distances between random sequences:
  ijs <- expand_grid(i = 1:m, j = 1:m) |>
    filter(i < j) |>
    slice_sample(n = n_rnd)
  
  results_rnd <- future_lapply(1:n_rnd, function(k) {
    hamming_rnd <- rep(NA, l)
    grantham_rnd <- rep(NA, l)
    try( {
      # sequence 1:
      i <- ijs$i[k]
      load(paste0(alig_path, "/", names(target_sequences)[i], ".RData"))
      coords1 <- coords_alig$coordinates
      rpoB_target1_AA <- translate(target_sequences[[i]])
      
      # sequence 2:
      j <- ijs$j[k]
      load(paste0(alig_path, "/", names(target_sequences)[j], ".RData"))
      coords2 <- coords_alig$coordinates
      rpoB_target2_AA <- translate(target_sequences[[j]])
      
      for(pos_Ecoli in 1:(l-1)) { # loop through E. coli rpoB
        pos_target1 <- translateCoordinate(pos_Ecoli, coords1, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
        pos_target2 <- translateCoordinate(pos_Ecoli, coords2, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
        if ((!is.na(pos_target1)) && (!is.na(pos_target2))) {
          AA_target1 <- as.character(rpoB_target1_AA[pos_target1])
          AA_target2 <- as.character(rpoB_target2_AA[pos_target2])
          hamming_rnd[pos_Ecoli] <- (AA_target1 != AA_target2)
          if ((AA_target1 %in% colnames(grantham)) && (AA_target2 %in% colnames(grantham))) {
            grantham_rnd[pos_Ecoli] <- grantham[AA_target1, AA_target2]
          } else {
            warning(paste0("Unknown AA (", AA_target, ")."))
          }
        } else if (xor(is.na(pos_target1), is.na(pos_target2))) {
          hamming_rnd[pos_Ecoli] <- TRUE
          grantham_rnd[pos_Ecoli] <- NA
        } else {
          hamming_rnd[pos_Ecoli] <- NA
          grantham_rnd[pos_Ecoli] <- NA
        }
      }
    }
    )
    return(list(hamming_rnd = hamming_rnd, grantham_rnd = grantham_rnd))
  }
  )
  
  hamming_Ecoli <- matrix(NA, nrow = m, ncol = l,
                          dimnames = list(names(target_sequences), paste0("P", 1:l)))
  grantham_Ecoli <- hamming_Ecoli
  for(i in 1:m) {
    hamming_Ecoli[i,] <- results_Ecoli[[i]]$hamming_Ecoli
    grantham_Ecoli[i,] <- results_Ecoli[[i]]$grantham_Ecoli
  }
  
  hamming_rnd <- matrix(NA, nrow = n_rnd, ncol = l,
                        dimnames = list(paste0(names(target_sequences)[ijs$i], 
                                               ":::",
                                               names(target_sequences)[ijs$j]),
                                        paste0("P", 1:l)))
  grantham_rnd <- hamming_rnd
  for(i in 1:n_rnd) {
    hamming_rnd[i,] <- results_rnd[[i]]$hamming_rnd
    grantham_rnd[i,] <- results_rnd[[i]]$grantham_rnd
  }
  
  # means aross all sequences:
  means = tibble(hamming_Ecoli = colMeans(hamming_Ecoli, na.rm = TRUE),
                 grantham_Ecoli = colMeans(grantham_Ecoli, na.rm = TRUE),
                 hamming_rnd = colMeans(hamming_rnd, na.rm = TRUE),
                 grantham_rnd = colMeans(grantham_rnd, na.rm = TRUE))
  
  return(list(hamming_Ecoli = hamming_Ecoli, 
              grantham_Ecoli = grantham_Ecoli,
              hamming_rnd = hamming_rnd, 
              grantham_rnd = grantham_rnd,
              means = means))
}



# Functions to generate "mutations tables", i.e. tables listing mutations
# in one or several sequences. Each row in such a table represents a particular
# mutation, and the fillMutationsTable function tries to make sure that the table
# is comprehensively filled.
#' Checks if a row in a mutations table indicates a point mutation
#'
#' @param muts A mutations table
#' @param k Row within that table.
#'
#' @return TRUE or FALSE
#'
#' @keywords internal
#'
isPointMutationNt <- function(muts, k) {
  pointMut <- TRUE

  # check nucleotides:
  if (!is.na(muts$Nt_OtoM[k])) {
    if (!stringr::str_detect(muts$Nt_OtoM[k], " -> "))
      pointMut <- FALSE
    nts <- stringr::str_split(muts$Nt_OtoM[k], pattern = stringr::fixed(" -> "))[[1]]
    if ((length(nts) != 2L) || (nchar(nts[1]) != 1L) || (nchar(nts[2]) != 1L))
      pointMut <- FALSE
  }

  if ((!is.na(muts$Nt_original[k])) && (nchar(muts$Nt_original[k]) != 1L)) {
    pointMut <- FALSE
  }

  if ((!is.na(muts$Nt_mutation[k])) && (nchar(muts$Nt_mutation[k]) != 1L)) {
    pointMut <- FALSE
  }

  return(pointMut)
}

# internal function to fill the Nt columns in a mutations table row 
# including Nt_original, Nt_mutation, Nt_OtoM, Nt_pos and Nt_pos_Ecoli:
fillMutationsTableRow_Nt <- function(muts, k, refSeq_ID, seqs, coordinates) {
  for(i in 1:2) { # repeat filling to make sure all columns are captured
    # fill Nt_OtoM from Nt_original + Nt_mutation:
    if (!is.na(muts$Nt_original[k]) && !is.na(muts$Nt_mutation[k])) {
      Nt_OtoM <- paste0(muts$Nt_original[k], " -> ", muts$Nt_mutation[k])
      if (is.na(muts$Nt_OtoM[k])) {
        muts$Nt_OtoM[k] <- Nt_OtoM
      } else {
        if (muts$Nt_OtoM[k] != Nt_OtoM) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_OtoM inconsistent with Nt_original and Nt_mutation"
        }
      }
    }

    # fill Nt_original and Nt_mutation from Nt_OtoM:
    if (!is.na(muts$Nt_OtoM[k])) {
      Nt_original <- stringr::word(muts$Nt_OtoM[k], 1)
      Nt_mutation <- stringr::word(muts$Nt_OtoM[k], 2, sep = stringr::fixed(" -> "))
      if (is.na(muts$Nt_original[k])) {
        muts$Nt_original[k] <- Nt_original
      } else {
        if (muts$Nt_original[k] != Nt_original) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_original inconsistent with Nt_OtoM"
        }
      }
      if (is.na(muts$Nt_mutation[k])) {
        muts$Nt_mutation[k] <- Nt_mutation
      } else {
        if (muts$Nt_mutation[k] != Nt_mutation) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_mutation inconsistent with Nt_OtoM"
        }
      }
    }

    # fill Nt_pos from Nt_pos_Ecoli:
    if (!is.na(muts$Nt_pos[k])) {
      Nt_pos_Ecoli <- translateCoordinate(muts$Nt_pos[k],
                                          coordinates[[refSeq_ID]],
                                          direction = "FocalToRef")
      if (is.na(Nt_pos_Ecoli)) {
        if (is.na(muts$Warning[k]))
          muts$Warning[k] <- "Nt_pos predicted to not exist in E. coli"
      }
      else if (is.na(muts$Nt_pos_Ecoli[k])) {
        muts$Nt_pos_Ecoli[k] <- Nt_pos_Ecoli
      } else {
        if (muts$Nt_pos_Ecoli[k] != Nt_pos_Ecoli) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_pos_Ecoli inconsistent with Nt_pos"
        }
      }
    }

    # fill Nt_pos_Ecoli from Nt_pos:
    if (!is.na(muts$Nt_pos_Ecoli[k])) {
      Nt_pos <- translateCoordinate(muts$Nt_pos_Ecoli[k],
                                    coordinates[[refSeq_ID]],
                                    direction = "RefToFocal")
      if (is.na(Nt_pos)) {
        if (is.na(muts$Warning[k]))
          muts$Warning[k] <- "Nt_pos_Ecoli predicted to not exist in this species"
      }
      else if (is.na(muts$Nt_pos[k])) {
        muts$Nt_pos[k] <- Nt_pos
      } else {
        if (muts$Nt_pos[k] != Nt_pos) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_pos inconsistent with Nt_pos_Ecoli"
        }
      }
    }

    # filling in Gene, Nt_original, Nt_pos and Nt_mutation from mut_name:
    if (!is.na(muts$Nt_mut_name[k])) {
      gene <- stringr::str_split(muts$Nt_mut_name[k], "_")[[1]][1]
      mutName <- stringr::str_split(muts$Nt_mut_name[k], "_")[[1]][2]
      Nt_pos <- as.integer(stringr::str_extract(mutName, "\\d+"))
      Nt_original <- stringr::str_extract_all(mutName, "[A-Z]+")[[1]][1]
      Nt_mutation <- stringr::str_extract_all(mutName, "[A-Z]+")[[1]][2]
      if (is.na(muts$Gene[k])) {
        muts$Gene[k] <- gene
      } else {
        if(muts$Gene[k] != gene) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_pos inconsistent with Nt_pos_Ecoli"
        }
      }
      if (is.na(muts$Nt_pos[k])) {
        muts$Nt_pos[k] <- Nt_pos
      } else {
        if(muts$Nt_pos[k] != Nt_pos) {
           if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_pos inconsistent with Nt_mut_name"
        }

      }
      if (is.na(muts$Nt_original[k])) {
        muts$Nt_original[k] <- Nt_original
      } else {
        if(muts$Nt_original[k] != Nt_original) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_pos inconsistent with Nt_mut_name"
        }
      }
      if (is.na(muts$Nt_mutation[k])) {
        muts$Nt_mutation[k] <- Nt_mutation
      } else {
        if (muts$Nt_mutation[k] != Nt_mutation) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_pos inconsistent with Nt_mut_name"
        }
      }
    }

    # filling Nt_mut_name from Nt_pos, Nt_original & Nt_mutation:
    if (!is.na(muts$Nt_pos[k]) && !is.na(muts$Nt_original[k]) && !is.na(muts$Nt_mutation[k])) {
      mutName <- paste0(muts$Gene[k], "_",
                        muts$Nt_original[k],
                        muts$Nt_pos[k],
                        muts$Nt_mutation[k])
      if (is.na(muts$Nt_mut_name[k])) {
        muts$Nt_mut_name[k] <- mutName
      } else {
        if (muts$Nt_mut_name[k] != mutName) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_mut_name inconsistent with Nt_original, Nt_pos and Nt_mutation"
        }
      }
    }

    # filling in Nt_original, Nt_pos_Ecoli and Nt_mutation from mut_name_Ecoli:
    if (!is.na(muts$Nt_mut_name_Ecoli[k])) {
      gene <- stringr::str_split(muts$Nt_mut_name_Ecoli[k], "_")[[1]][1]
      mutName_Ecoli <- stringr::str_split(muts$Nt_mut_name_Ecoli[k], "_")[[1]][2]
      Nt_pos_Ecoli <- as.integer(stringr::str_extract(mutName_Ecoli, "\\d+"))
      Nt_original <- stringr::str_extract_all(mutName_Ecoli, "[A-Z]+")[[1]][1]
      Nt_mutation <- stringr::str_extract_all(mutName_Ecoli, "[A-Z]+")[[1]][2]
      if (is.na(muts$Nt_pos_Ecoli[k])) {
        muts$Nt_pos_Ecoli[k] <- Nt_pos_Ecoli
      } else {
        if(muts$Nt_pos_Ecoli[k] != Nt_pos_Ecoli) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_pos_Ecoli inconsistent with Nt_mut_name_Ecoli"
        }
      }
      if (is.na(muts$Nt_original[k])) {
        muts$Nt_original[k] <- Nt_original
      } else {
        if(muts$Nt_original[k] != Nt_original) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_original inconsistent with Nt_mut_name_Ecoli"
        }
      }
      if (is.na(muts$Nt_mutation[k])) {
        muts$Nt_mutation[k] <- Nt_mutation
      } else {
        if(muts$Nt_mutation[k] != Nt_mutation) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_mutation inconsistent with Nt_mut_name_Ecoli"
        }
      }
    }

    # filling Nt_mut_name_Ecoli from Nt_pos, Nt_original & Nt_mutation:
    if (!is.na(muts$Nt_pos_Ecoli[k]) && !is.na(muts$Nt_original[k]) && !is.na(muts$Nt_mutation[k])) {
      mutName_Ecoli <- paste0(muts$Gene[k], "_",
                              muts$Nt_original[k],
                              muts$Nt_pos_Ecoli[k],
                              muts$Nt_mutation[k])
      if (is.na(muts$Nt_mut_name_Ecoli[k])) {
        muts$Nt_mut_name_Ecoli[k] <- mutName_Ecoli
      } else {
        if (muts$Nt_mut_name_Ecoli[k] != mutName_Ecoli) {
          if (is.na(muts$Warning[k]))
            muts$Warning[k] <- "Nt_mut_name_Ecoli inconsistent with Nt_original, Nt_pos_Ecoli and Nt_mutation"
        }
      }
    }
  }
  return(muts)
}

#' Fill in missing data in a mutations table.
#'
#' This function tries to infer missing data in a mutations table,
#' focusing on nucleotide-level information. It resolves explicit positions
#' of mutations where this data is only implicitly available and standardizes
#' positions to reference (E.coli) coordinates.
#' The result is a standardised table that can then more readily be analysed.
#'
#' @param muts An (incomplete) mutations table.
#' @param refs A table specifying the species and names for all reference sequences.
#' @param seqs Reference sequences as a DNAStringSet object, with names as specified in `refs`.
#' @param coordinates A list of coordinates, with names corresponding to seqs.
#'
#' @return A muts table that is (hopefully) more complete than the input table.
#' @export
#'
fillMutationsTableNt <- function(muts, refs, seqs, coordinates) {
  # add additional columns:
  muts <- dplyr::mutate(muts, Strain_ID = gsub(" ", "_", paste(Species, Strain)),
                              RefSeq_ID = NA,
                              Warning = NA)

  # check if pos columns can be converted to numbers, if not flag warning:
  muts$Warning[grepl("\\D", muts$Nt_pos)] <- "Nt_pos not a number"
  muts$Warning[grepl("\\D", muts$Nt_pos_Ecoli)] <- "Nt_pos_Ecoli not a number"
  # muts$Warning[grepl("\\D", muts$AA_pos)] <- "AA_pos not a number"
  # muts$Warning[grepl("\\D", muts$AA_pos_Ecoli)] <- "AA_pos_Ecoli not a number"

  # convert positions into numbers:
  # (needs to be improved if mutations other than point mutations are to be included)
  muts$Nt_pos <- suppressWarnings(as.integer(muts$Nt_pos))
  muts$Nt_pos_Ecoli <- suppressWarnings(as.integer(muts$Nt_pos_Ecoli))
  # muts$AA_pos <- suppressWarnings(as.integer(muts$AA_pos))
  # muts$AA_pos_Ecoli <- suppressWarnings(as.integer(muts$AA_pos_Ecoli))
  #muts$MIC_mgPerL <- as.double(muts$MIC_mgPerL)

  # infer gene name from mutation name:
  # for(i in 1:nrow(muts)) {
  #   if (is.na(muts$Gene[i])) {
  #     mut_names <- c(muts$Nt_mut_name[i], muts$Nt_mut_name_Ecoli[i],
  #                    muts$AA_mut_name[i], muts$AA_mut_name_Ecoli[i])
  #     mut_names <- mut_names[!is.na(mut_names)]
  #     if (length(mut_names) > 0) {
  #       gene_names <- str_extract(mut_names, "^.+?(?=_)")
  #       if (length(unique(gene_names)) > 1) {
  #         muts$Warning[i] <- "inconsistent gene names"
  #       } else {
  #         muts$Gene[i] <- gene_names[1]
  #       }
  #     } else {
  #       muts$Warning[i] <- "missing gene name"
  #     }
  #   }
  # }

  genes <- unique(muts$Gene)
  strainIDs <- unique(muts$Strain_ID)
  for(i in 1:length(genes)) {
    for(j in 1:length(strainIDs)) {
      refSeq_ID <- NA
      # check if there are any mutations in that strain j and gene i:
      if ((muts |>
           dplyr::filter(Gene == genes[i], Strain_ID == strainIDs[j]) |>
           nrow()) > 0) {
        geneStrain <- paste(genes[i], strainIDs[j], sep = "_")
        cat(paste0("Working on gene_species_strain ", geneStrain, ".\n"))
        # finding the appropriate reference sequence for the gene:
        if (geneStrain %in% names(seqs)) {  # specific strain has reference
          refSeq_ID <- geneStrain
          cat("  -> Specific reference for this strain found.\n")
        } else {
          species <- paste(stringr::str_split(strainIDs[j], "_")[[1]][1:2],
                           collapse = " ")
          filteredRefs <- dplyr::filter(refs,
                                        Gene == genes[i] & Species == species & RefStrain)
          if ((nrow(filteredRefs) == 1) && (filteredRefs$FASTA_name[1] %in% names(seqs))) {
            refSeq_ID <- filteredRefs$FASTA_name[1]
            cat(paste0("  -> No specific reference found, but used ", filteredRefs$FASTA_name[1], " instead.\n"))
          } else {
            cat("  -> No reference sequence found.\n")
          }
        }
      }
      if (!is.na(refSeq_ID)) {
        for(k in 1:nrow(muts)) {
          if ((muts$Gene[k] == genes[i]) &&
              (muts$Strain_ID[k] == strainIDs[j])) {
            muts$RefSeq_ID[k] <- refSeq_ID
            if (all(is.na(c(muts$Nt_pos[k], muts$Nt_pos_Ecoli[k],
                            muts$Nt_mut_name[k], muts$Nt_mut_name_Ecoli[k],
                            muts$AA_pos[k], muts$AA_pos_Ecoli[k],
                            muts$AA_mut_name[k], muts$AA_mut_name_Ecoli[k]))) && (is.na(muts$Warning[k]))) {
              muts$Warning[k] <- "No mutation position (Nt or AA) available"
            }
            if (all(is.na(c(muts$Nt_mutation[k],
                            muts$Nt_mut_name[k], muts$Nt_mut_name_Ecoli[k],
                            muts$Codon_mutation[k],
                            muts$AA_mutation[k],
                            muts$AA_mut_name[k], muts$AA_mut_name_Ecoli[k]))) && (is.na(muts$Warning[k]))) {
              muts$Warning[k] <- "No mutation (Nt, codon or AA) available"
            }
            if (is.na(muts$Warning[k])) {
              if (isPointMutationNt(muts, k)) {
                muts <- fillMutationsTableRow_Nt(muts, k, refSeq_ID, seqs, coordinates)
              } else {
                if (is.na(muts$Warning[k]))
                  muts$Warning[k] <- "Not a single point mutation"
              }
            }
          }
        }
      }
    }
  }
  muts$Warning[is.na(muts$RefSeq_ID)] <- "No reference sequence available"
  mutationsTableSummaryNt(muts)
  return(muts)
}


#' Print a summary of a mutations table
#'
#' @param muts A mutations table, as created using the fillMutationsTableNt function.
#'
#' @export
#'
mutationsTableSummaryNt <- function(muts) {
  cat("\nSummary of mutations table:\n")
  cat("  Number of species:",
      muts |> dplyr::pull(Species) |> unique() |> length(), "\n")
  cat("  Genes with mutations:",
      paste(muts |>
              dplyr::pull(Gene) |>
              unique(),
            collapse = ", "),
      "\n")
  cat("  Total number of recorded mutations:", nrow(muts), "\n")
  cat("  Number of unique Nt mutations across species:",
      muts |>
        dplyr::filter(!is.na(Nt_mut_name)) |>
        dplyr::group_by(Species, Nt_mut_name) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
        dplyr::pull(n) |>
        sum(),
      "\n")
  cat("  Number of unique orthologous Nt mutations:",
      muts |>
        dplyr::filter(!is.na(Nt_mut_name)) |>
        dplyr::pull(Nt_mut_name_Ecoli) |>
        unique() |>
        length(),
      "\n")

  missingRefs <- muts |>
    dplyr::filter(is.na(RefSeq_ID)) |>
    dplyr::select(Gene, Strain_ID) |>
    dplyr::distinct(Gene, Strain_ID) |>
    tidyr::unite(Gene_Strain, Gene, Strain_ID, sep = "_") |>
    dplyr::pull(Gene_Strain)
  if (length(missingRefs) > 0) {
    cat("\nNo reference found for the following Gene_Species_Strain combinations:\n")
    cat(paste(" ", missingRefs, sep = "", collapse = "\n"))
    cat("\n")
  }

  nWarnings <- sum(!is.na(muts$Warning))
  if (nWarnings == 1L) {
    cat(paste0("\nProblems encountered in one row of the table (see Warnings column).\n\n"))
  } else if (nWarnings > 0) {
    cat(paste0("\nProblems encountered in ", nWarnings, " rows of the table (see Warnings column).\n\n"))
  }
}


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
  coords_alig <- ALJEbinf::getCoordinates(target_sequence, reference_Ecoli, aligOutput = TRUE)
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

#' filters final_output to retrieve data related to high quality sequences and reliable mutation list
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


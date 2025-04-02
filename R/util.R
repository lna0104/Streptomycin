# function to extract the species name from the target sequence name 
# (e.g., "rpoB_1_Escherichia_coli_GCF_000005845.2" --> "Escherichia coli")
extract_species_name <- function(target_name, target_gene) {
  #return(gsub("rpoB_([A-z0-9_()\\.\\'-]+)_GCF_.+", "\\1", target_name))
  return(gsub(paste0(target_gene, "_[0-9]+_([A-z0-9_ ()\\.\\'-]+)_GCF_.+"), "\\1", target_name))
}

# function to extract the accession number from the target sequence name 
# (e.g., "rpoB_1_Escherichia_coli_GCF_000005845.2" --> "GCF_000005845.2")
extract_accession_number <- function(target_name, target_gene) {
  return(gsub(paste0(target_gene, "_[0-9]+_[A-z0-9_ ()\\.\\'-]+_(GCF_.+)"), "\\1", target_name))
}

# function to extract the accession number from the target sequence name 
# (e.g., "rpoB_1_Escherichia_coli_GCF_000005845.2" --> 1)
extract_gene_copy_number <- function(target_name, target_gene) {
  return(gsub(paste0(target_gene, "_([0-9]+)[A-z0-9_ ()\\.\\'-]+_GCF_.+"), "\\1", target_name))
}


#' Compute the distance between two sequences 
#'
#' @param seqs an object of class DNAStringSet as the subject sequence (target)
#' @param ref an object of class DNAStringSet as the pattern sequence (reference)
#'
#' @return a data frame providing the distance of the subject against the pattern
#' @export
#'
#' @examples get_dist(rpoB_target_sequences, rpoB_reference_E.coli)
get_dist <- function(seqs, ref){
  distance_check <- data.frame(seq_name = extract_species_name(names(seqs)), 
                               dist_from_ref = NA)
  for(i in 1:nrow(distance_check)) {
    dis <- Biostrings::stringDist(Biostrings::translate(c(seqs[i], ref), if.fuzzy.codon = "X"),
                                  method = "levenshtein")
    distance_check$dist_from_ref[i] <- dis
    Sys.sleep(0.1)
    cat("\rScreened sequences:", i, "/", length(seqs))
    flush.console() 
    cat("\n")
  }
  return(distance_check)
}


#' retrieve three-letter nucleotide codon of an amino acid in a AAString based on corresponding DNA sequence
#'
#' @param dna_seq a DNAString object 
#' @param AA_pos a numeric vector providing the position of amino acid
#'
#' @return a character string showing three-letter nucleotide codon of an amino acid in a AAString
#' @export
#'
#' @examples get_nucleotide_codon(rpoB_seq, 501)
get_nucleotide_codon <- function(dna_seq, AA_pos){
  AA_pos <- AA_pos
  nt_codon <- as.character(substr(dna_seq, 
                                  AA_pos*3 - 2, 
                                  AA_pos*3))
  return(nt_codon)
}



# function to generate alignment scores of random sequences:
rAlignmentScore <- function(n, seq, l = length(seq)) {
  scores <- rep(NA, n)
  for(i in 1:n) {
    scores[i] <- sample(c("A", "T", "C", "G"), l, replace = TRUE) |>
      paste0(collapse = "") |>
      DNAString() |>
      Biostrings::translate() |>
      pairwiseAlignment(Biostrings::translate(seq)) |>
      Biostrings::score()
  }
  return(scores)
}

# example:
# seq <- readDNAStringSet("./data/rpoB_references.fasta")[["rpoB_Escherichia_coli_MG1655"]]
# rAlignmentScore(10, seq, l = 600)

get_scores <- function(seqs, ref) {
  align <- lapply(seqs, function(seq) pairwiseAlignment(ref, seq))
  scores <- sapply(align, function(aln) aln@score)
  info <- data.frame(scores = scores, species = names(seqs))
  return(info)
}

# Hamming distance between two vectors:
hamming <- function(v1, v2) {
  return(sum(v1 != v2))
}

# Hamming distance between two strings of equal length:
# The function is parallelised so that x and y can also be vectors of strings.
# (Recycling is not supported though at this point.)
hamming_str <- function(x, y) {
  x_split <- str_split(x, pattern = '')
  y_split <- str_split(y, pattern = '')
  return(sapply(1:length(x), function(i) {hamming(x_split[[i]], y_split[[i]])} ))
}


#' Run Expression, Retrying On Error
#'
#' Runs a given expression repeatedly.
#' If it errors, prints the error and reruns the expression.
#' Useful for web requests that might fail intermittently.
#'
#' @param expr an R expression to try.
#' @param max_tries the maximum number of times to try the expression.
#' @param envir the environment in which to run the expression.
#'
#' @returns the return value of the expression on its successful execution
#' @export
#'
#' @examples
#' retry(stopifnot(runif(1) > 0.95))
retry <- function(expr, max_tries = Inf, envir = parent.frame()) {
  expr <- substitute(expr)
  while(max_tries > 1) {
    result <- tryCatch(eval(expr, envir), error = function(e) {return(e)})
    if(!is(result, "error")) {
      return(result)
    } else {
      cat(sprintf("Error in %s : %s\n", deparse(result$call), result$message), file = stderr())
      cat("Retrying...\n", file = stderr())
      max_tries <- max_tries - 1
    }
  }
  #If this is our last try, just run without error handling.
  eval(expr, envir)
}


#' Find papers where particular mutations have been reported in
#'
#' @param mut_names A vector of one or more mutation names (e.g., "526C")
#' @param muts Completed table of reported mutations
#' @param reports_references Table of bibiographical information
#'
#' @return A table with all papers for each of the querried mutation.
#' @export
#'
find_papers <- function(mut_names, muts, reports_references) {
  papers <- muts |>
    mutate(mut_name = paste0(AA_pos_Ecoli, AA_mutation)) |>
    filter(mut_name %in% mut_names) |>
    left_join(reports_references, by = join_by(Species, Ref_code)) |>
    select(mut_name, Species, Origin, Ref_code, Full_ref, DOI) |>
    distinct()
  return(papers)
}


# Split interval 1:n into m roughly equal intervals:
# (courtesy of chatGPT)
split_indices <- function(n, m) {
  groups <- cut(seq_len(n), breaks = m, labels = FALSE)
  data.frame(
    start = tapply(seq_len(n), groups, FUN = function(x) x[1]),
    end = tapply(seq_len(n), groups, FUN = function(x) x[length(x)])
  )
}
#' Tests if tips of a tree are more closely to each other than expected by chance
#'
#' @param tree A tree is phylo format
#' @param tips A vector of tip names
#' @param n Number of permutations
#'
#' @return A p values, giving the probability that the tree has the
#' same or a smaller mean distance between permutated tips than
#' with the actually tested tips.
#' 
relatedness_test <- function(tree, tips, n = 1000) {
  d_matrix <- cophenetic.phylo(tree)
  d_true <- mean(d_matrix[tips, tips])
  d_permut <- replicate(n, { 
    rnd_tips <- sample(tree$tip.label, length(tips))
    mean(d_matrix[rnd_tips, rnd_tips])
  })
  # p value, calculated as percentile of true value within permutated d values:
  p_value <- ecdf(d_permut)(d_true)
  result <- list(mean_distance_tips = d_true,
                 mean_distance_tips_permutated = d_permut,
                 P = p_value)
  return(result)
}


get_phylosignals_nt <- function(subtree, 
                             species_output,
                             sample_n = NULL) {
  ##change tip labels to match them to data:  
  subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))
  
  # subsampling from tree, meant for testing only:
  if (!is.null(sample_n)) {
    m <- length(subtree$tip.label) - sample_n
    subtree <- drop.tip(subtree, sample(subtree$tip.label)[1:m])
  }
  
  # retain only species that are in the tree:
  species_output <- species_output |>
    filter(species %in% subtree$tip.label)
  
  # resistant species:
  resistant_species <- species_output |>
    filter(resistance == "resistant") |>
    pull(species)
  
  # obtain evolvabilities as named vectors:
  evolvabilityI <- species_output |>
    select(species, evolvabilityI) |>
    deframe()
  # evolvabilityII <- species_output |>
  #   select(species, evolvabilityII) |>
  #   deframe()
  
  cat("Performing permutation test for intrinsic resistance...")
  permtest_resistance <- relatedness_test(subtree, resistant_species)
  # cat("done!\nCalculating Blomberg's K for evolvability I...")
  # K_evolvabilityI <- phytools::phylosig(subtree, evolvabilityI, method = "K", test = TRUE)
  # cat("done!\nCalculating Pagel's lambda for evolvability I...")
  # lambda_evolvabilityI <- phytools::phylosig(subtree, evolvabilityI, method = "lambda", test = TRUE)
  # cat("done!\nCalculating Blomberg's K for evolvability II...")
  # K_evolvabilityII <- phytools::phylosig(subtree, evolvabilityII, method = "K", test = TRUE)
  # cat("done!\nCalculating Pagel's lambda for evolvability II.\n")
  # lambda_evolvabilityII <- phytools::phylosig(subtree, evolvabilityII, method = "lambda", test = TRUE)
  
  return(list(permtest_resistance = permtest_resistance
              # K_evolvabilityI = K_evolvabilityI,
              # K_evolvabilityII = K_evolvabilityII,
              # lambda_evolvabilityI = lambda_evolvabilityI,
              # lambda_evolvabilityII = lambda_evolvabilityII
              ))
}
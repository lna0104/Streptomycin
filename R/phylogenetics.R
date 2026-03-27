#' build a tree based on GTDB bacterial life tree using species names and accessions
#'
#' @param output a data frame containing the predicted mutations for all species
#' @param original_tree a phylo object providing the bacterial tree of life (from GTDB)
#' @param meta_data a data frame providing the information of species which are included in bacterial tree of life (from GTDB)
#'
#' @return a list providing the subset tree and a data frame presenting all tip labels
#' @export
#'
#' @examples get_subtree(filtered_output, original_tree, meta_data, outliers)
#' #
get_subtree <- function(output, original_tree, meta_data, outliers) {
  #  our_species <- output |>
  #    select(species, accession_numbers) |>
  #    distinct()
  #  #1. create a column "species_accessions" in metadata
  #  meta_data_mod <- meta_data %>%
  #   select(accession, ncbi_genbank_assembly_accession, ncbi_organism_name) %>%
  #   mutate(species_name = gsub("\\[|\\]", "", ncbi_organism_name)) %>%
  #   mutate(species_name = gsub("\\'|\\'", "", species_name)) %>%
  #   mutate(species_name = str_replace_all(species_name, "Candidatus", "")) %>%
  #   mutate(species_name = gsub("^ ", "", species_name)) %>%
  #   mutate(species_name = sub("^([^ ]+[ ]+[^ ]+).*$", "\\1", #keeps only the first two words of a character element
  #                             species_name)) %>%
  #   filter(accession %in% original_tree$tip.label) %>%
  #   mutate(acc_ID = paste(ncbi_genbank_assembly_accession,
  #                         species_name, sep = "_"))
  #  #2. change tip labels from "accessions" to "species_accessions"
  #  original_tree$tip.label <- unname(setNames(meta_data_mod$acc_ID, meta_data_mod$accession)[original_tree$tip.label])
  #  original_tree$tip.label <- paste0(sub("A", "F", original_tree$tip.label))
  #  #3. exploring our species in metadata
  #  accession_not_tree <- setdiff(our_species$accession_numbers, substr(original_tree$tip.label, -1, 15))
  #  species_accession_not_tree <- our_species %>%
  #    filter(accession_numbers %in% accession_not_tree) %>%
  #    pull(species)
  #  #4. subset original tree based on tip labels
  #  #subset the original tree based on presence of species names or accessions in tips:
  #  subset_tip_indices_species_name <- sapply(species_accession_not_tree, function(subset_label) {
  #    which(grepl(subset_label, original_tree$tip.label, fixed = TRUE))[1] #checks the presence of each species names in all tip labels, and returns index of first one
  #  })
  #  #subset based on accessions
  #  subset_tip_indices_accession <- sapply(our_species$accession_numbers, function(acc) {
  #    which(grepl(acc, original_tree$tip.label, fixed = TRUE))[1] #checks the presence of each species names in all tip labels, and returns index of first one
  #  })
  #  #subset for missing species
  #  # subset_tip_indices_missing_names <- sapply(missing_names$names_original_tree, function(subset_label) {
  #  #   which(grepl(subset_label, original_tree$tip.label, fixed = TRUE))[1] #checks the presence of each species names in all tip labels, and returns index of first one
  #  # })
  #  #subset for outliers to be removed
  #  subset_tip_indices_outliers <- which(purrr::reduce(lapply(outliers$outliers, function(subset_label) {
  #    grepl(subset_label, original_tree$tip.label, fixed = TRUE) #checks the presence of each species names in all tip labels
  #  }), `|`))
  #  #sum of all subsets to be kept
  #  subset_tip_indices <- unique(c(subset_tip_indices_species_name,
  #                                 subset_tip_indices_accession))
  #  subset_tip_indices <- subset_tip_indices[!is.na(subset_tip_indices)]
  #  #6. get the sub_tree
  #  subset_tree <- get_subtree_with_tips(original_tree, only_tips = subset_tip_indices,
  #                                       omit_tips = subset_tip_indices_outliers, force_keep_root = TRUE)$subtree
  #  ### I did some extra coding to find actual name of those tip labels which only have "genus name sp."
  #  #make a data frame using subset_tree$tip.label
  #  tip_labels <- data.frame(old_tip_label = subset_tree$tip.label) %>%
  #    mutate(species_phylo = sub(".*_", "", old_tip_label)) %>%
  #    mutate(acc = substr(old_tip_label, 1, 15)) %>% #separate accession numbers as a column
  #    left_join(our_species, by = join_by(acc == accession_numbers)) %>%
  #    select(acc, species_phylo, species) %>%
  #    distinct() %>%
  #    mutate(species = ifelse(is.na(species), species_phylo, species)) %>%
  #    mutate(new_tip_label = paste0(acc, "_", species))
  #  #replace new names in tip.labels of tree
  #  subset_tree$tip.label <- tip_labels$new_tip_label
  #  return(list(tree = subset_tree,
  #              tips = tip_labels))

  our_species <- output |>
    select(species, accession_numbers) |>
    distinct()

  # 1. create a column "species_accessions" in metadata
  meta_data_tree <- meta_data %>%
    select(accession, ncbi_genbank_assembly_accession, species) %>%
    filter(accession %in% original_tree$tip.label) %>%
    mutate(acc_ID = paste(ncbi_genbank_assembly_accession,
      species,
      sep = "_"
    ))

  # 2. change tip labels from "accessions" to "species_accessions"
  original_tree$tip.label <- unname(setNames(meta_data_tree$acc_ID, meta_data_tree$accession)[original_tree$tip.label])
  original_tree$tip.label <- paste0(sub("A", "F", original_tree$tip.label))

  # 3. exploring our species in metadata

  accession_not_tree <- setdiff(our_species$accession_numbers, substr(original_tree$tip.label, -1, 15))

  species_accession_not_tree <- our_species %>%
    filter(accession_numbers %in% accession_not_tree) %>%
    pull(species)

  # 4. subset original tree based on tip labels

  # subset the original tree based on presence of species names or accessions in tips:
  subset_tip_indices_species_name <- sapply(species_accession_not_tree, function(subset_label) {
    which(grepl(subset_label, original_tree$tip.label, fixed = TRUE))[1] # checks the presence of each species names in all tip labels, and returns index of first one
  })

  # subset based on accessions
  subset_tip_indices_accession <- sapply(our_species$accession_numbers, function(acc) {
    which(grepl(acc, original_tree$tip.label, fixed = TRUE))[1] # checks the presence of each species names in all tip labels, and returns index of first one
  })

  # subset for outliers to be removed
  # subset_tip_indices_outliers <- which(purrr::reduce(lapply(outliers$outliers, function(subset_label) {
  #   grepl(subset_label, original_tree$tip.label, fixed = TRUE) # checks the presence of each species names in all tip labels
  # }), `|`))

  # sum of all subsets to be kept
  subset_tip_indices <- unique(c(
    subset_tip_indices_species_name,
    subset_tip_indices_accession
  ))
  subset_tip_indices <- subset_tip_indices[!is.na(subset_tip_indices)]

  # 6. get the sub_tree
  subset_tree <- get_subtree_with_tips(original_tree,
    only_tips = subset_tip_indices,
    force_keep_root = TRUE
  )$subtree

  tip_labels <- data.frame(old_tip_label = subset_tree$tip.label) %>%
    mutate(
      acc = str_extract(old_tip_label, "^[^_]+_[^_]+"),
      species_phylo = str_remove(old_tip_label, "^[^_]+_[^_]+_")
    ) %>%
    left_join(our_species, by = join_by(acc == accession_numbers)) %>%
    select(acc, species_phylo, species) %>%
    distinct() %>%
    mutate(species = ifelse(is.na(species), species_phylo, species)) %>%
    mutate(new_tip_label = paste0(acc, "_", species))
  # replace new names in tip.labels of tree
  subset_tree$tip.label <- tip_labels$new_tip_label

  return(list(
    tree = subset_tree,
    tips = tip_labels
  ))
}


#' Generate a phylogenetic tree with bootstrapping values
#'
#' @param seqs a XStringset object
#' @return an object of class matrix which represents the binary matrix derived from the phylogenetic tree.
#' and at the same time returns the bootstrapping values for each node
#' @export
#'
#' @examples get_tree_with_boots(rpsL_target_sequences)
get_tree_with_boots <- function(seqs) {
  # keep long enough sequences
  seqs <- seqs[lengths(seqs) >= 2000]
  # perform multiple sequence alignment
  aligned_sequences <- muscle::muscle(seqs)
  # Convert alignment to a phyDat object
  phy_data <- as.phyDat(aligned_sequences)

  # compute all possible types of a given model for a given alignment
  estimated_models <- modelTest(phy_data,
    model = "K80",
    multicore = TRUE,
    mc.cores = detectCores()
  ) # count number of system cores
  # Generate or infer a tree using maximum likelihood (ML)
  computed_data <- pml_bb(estimated_models,
    rearrangement = "stochastic"
  )
  return(computed_data)
}


#' Tests if tips of a tree are more closely to each other than expected by chance
#'
#' @param tree A tree is phylo format
#' @param tips A vector of tip names
#' @param n Number of permutations
#'
#' @return A p values, giving the probability that the tree has the
#' same or a smaller mean distance between permutated tips than
#' with the actually tested tips.

relatedness_test <- function(tree, tips, n = 1000) {
  d_matrix <- cophenetic.phylo(tree)
  d_true <- mean(d_matrix[tips, tips])
  d_permut <- replicate(n, {
    rnd_tips <- sample(tree$tip.label, length(tips))
    mean(d_matrix[rnd_tips, rnd_tips])
  })
  # p value, calculated as percentile of true value within permutated d values:
  p_value <- ecdf(d_permut)(d_true)
  result <- list(
    mean_distance_tips = d_true,
    mean_distance_tips_permutated = d_permut,
    P = p_value
  )
  return(result)
}


get_phylosignals <- function(subtree,
                             species_output,
                             sample_n = NULL) {
  ## change tip labels to match them to data:
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
  evolvabilityII <- species_output |>
    select(species, evolvabilityII) |>
    deframe()

  cat("Performing permutation test for intrinsic resistance...")
  permtest_resistance <- relatedness_test(subtree, resistant_species)
  cat("done!\nCalculating Blomberg's K for evolvability I...")
  K_evolvabilityI <- phytools::phylosig(subtree, evolvabilityI, method = "K", test = TRUE)
  cat("done!\nCalculating Pagel's lambda for evolvability I...")
  lambda_evolvabilityI <- phytools::phylosig(subtree, evolvabilityI, method = "lambda", test = TRUE)
  cat("done!\nCalculating Blomberg's K for evolvability II...")
  K_evolvabilityII <- phytools::phylosig(subtree, evolvabilityII, method = "K", test = TRUE)
  cat("done!\nCalculating Pagel's lambda for evolvability II.\n")
  lambda_evolvabilityII <- phytools::phylosig(subtree, evolvabilityII, method = "lambda", test = TRUE)

  return(list(
    permtest_resistance = permtest_resistance,
    K_evolvabilityI = K_evolvabilityI,
    K_evolvabilityII = K_evolvabilityII,
    lambda_evolvabilityI = lambda_evolvabilityI,
    lambda_evolvabilityII = lambda_evolvabilityII
  ))
}


#' plot subtree
#'
#' @param subtree a phylo object of the subset tree
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpoB sequences across bacterial species
#' @param bacterial_taxonomy a data frame providing info on bacterial taxonomic phylogeny
#' @param outliers a data frame providing species which are misplaced inside the original tree
#' @param file_name the path that the plot should be save in
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#'
#' @export
#'
#' @examples plot_subtree(subtree, filtered_output, bacterial_taxonomy, "./plots/myfilename.pdf")
plot_subtree <- function(subtree, species_output, bacterial_taxonomy, file_name) {
  # clades to be labeled in the tree:
  clade_labels <- tibble::tribble(
    ~clade, ~tier,
    "Planctomycetia", "major",
    "Actinomycetes", "major",
    "Alphaproteobacteria", "major",
    "Sphingomonadales", "minor",
    "Devosiaceae", "minor",
    "Rickettsiales", "minor",
    "Coriobacteriia", "major",
    "Gammaproteobacteria", "major",
    "Bacteroidia", "major"
  )

  # change tip labels to make manipulations easier:
  subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))


  species_data <- species_output |>
    left_join(
      bacterial_taxonomy,
      by = join_by(genus),
      relationship = "many-to-many"
    ) |>
    mutate(
      major_clade = case_when(
        family %in% clade_labels$clade[clade_labels$tier == "major"] ~ family,
        order %in% clade_labels$clade[clade_labels$tier == "major"] ~ order,
        class %in% clade_labels$clade[clade_labels$tier == "major"] ~ class,
        phylum %in% clade_labels$clade[clade_labels$tier == "major"] ~ phylum,
        TRUE ~ "none"
      )
    )

  data_tree <- as_tibble(subtree) |>
    left_join(species_data, by = join_by(label == species)) |>
    filter(!is.na(node), !is.na(parent), !is.na(label)) |>
    replace_na(list(
      family = "undefined",
      order = "undefined",
      class = "undefined",
      major_clade = "none"
    ))

  ## identify common ancestor nodes for species in major bacterial orders
  get_nodes_for_clade <- function(dt, clade) {
    nodes <- c(
      dt$node[dt$phylum == clade],
      dt$node[dt$class == clade],
      dt$node[dt$order == clade],
      dt$node[dt$family == clade]
    )
    unique(na.omit(nodes))
  }

  clades_all <- clade_labels |>
    rowwise() |>
    mutate(
      nodes = list(get_nodes_for_clade(data_tree, clade)),
      common_ancestor = if (length(nodes) >= 2) ape::getMRCA(subtree, nodes) else NA_integer_
    ) |>
    ungroup() |>
    filter(!is.na(common_ancestor))

  major_df <- clades_all |> filter(tier == "major")
  minor_df <- clades_all |> filter(tier == "minor")

  ## make the final tree data as a S4 object for plotting
  subtree_data <- tidytree::as.treedata(data_tree)

  ## plotting the tree
  p_tree <- ggtree(subtree_data,
    layout = "circular",
    branch.length = "none",
    color = "#eee0cd"
  )

  # add predicted resistance:
  p1 <- gheatmap(p_tree, species_data |> select(species, resistance) |> column_to_rownames(var = "species"),
    offset = 0, color = NULL, colnames = FALSE, width = 0.05
  ) +
    scale_fill_manual(
      values = c("resistant" = "#ba181b", "susceptible" = "#d3d3d3"),
      labels = c("resistant", "susceptible", ""),
      na.value = hsv(0, 0, 0, 0),
      name = "Predicted\nresistance"
    )

  # add predicted evolvability:
  p2 <- gheatmap(p1 + new_scale_fill(),
    species_data |> select(species, evolvabilityI) |> column_to_rownames(var = "species"),
    offset = 3.4, color = NULL, colnames = FALSE, width = 0.05
  ) +
    scale_fill_gradientn(colours = c("#0077b6", "#ffd700"), na.value = "grey20", name = "Predicted\nevolvability") +
    scale_y_continuous(limits = c(0, round(1.005 * (nrow(data_tree) + 1) / 2)))

  # add major clade labels:
  p3 <- p2
  if (nrow(major_df) > 0) {
    major_cols <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Pastel1"))[1:nrow(major_df)]
    p3 <- p3 +
      geom_cladelab(
        node = major_df$common_ancestor,
        label = major_df$clade,
        align = TRUE,
        geom = "text",
        fontsize = 3.0,
        barcolour = major_cols,
        textcolour = major_cols,
        offset.text = 2.0,
        offset = 9,
        barsize = 2,
        horizontal = FALSE,
        hjust = 0.5
      )
  }

  # add minor clade labels
  if (nrow(minor_df) > 0) {
    minor_cols <- rep(brewer.pal(8, "Set2"), length.out = nrow(minor_df))
    for (i in seq_len(nrow(minor_df))) {
      p3 <- p3 +
        geom_cladelab(
          node = minor_df$common_ancestor[i],
          label = minor_df$clade[i],
          align = TRUE,
          geom = "text",
          fontsize = 3.0,
          barcolour = minor_cols[i],
          textcolour = minor_cols[i],
          offset.text = 2.0,
          offset = 11.5,
          barsize = 2,
          horizontal = FALSE,
          hjust = 0.5
        )
    }
  }

  p3 <- p3 + theme(legend.position = "bottom")

  ggsave(filename = file_name, p3, width = 10, height = 10)
}


#' plot subtree for a specific clade (genus, family, order, class)
#'
#' @param subtree a phylo object of the subset tree
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpsL
#'  sequences across bacterial species
#' @param bacterial_taxonomy a data frame providing info on bacterial taxonomic phylogeny
#' @param genus the name of the genus
#' @param family the name of the family
#' @param order the name of the order
#' @param class the name of the class
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#' @export
#'
#' @examples plot_subtree(subtree, filtered_output, bacterial_taxonomy, "./plots/myfilename.pdf")
plot_subtree_clade <- function(subtree,
                               species_output,
                               bacterial_taxonomy,
                               genus = NULL, family = NULL, order = NULL, class = NULL,
                               file_name,
                               tip_text_size = 4,
                               title_text_size = 18,
                               base_text_size = 14,
                               legend_text_size = 12,
                               tree_line_size = 0.4,
                               heatmap_width = 0.03,
                               fig_width = 16, # output width
                               tip_height = 0.25) { # height per tip label

  if (sum(!is.null(c(genus, family, order, class))) != 1L) {
    stop("One and only one of genus, family, order or class must be specified, the others need to be NULL.")
  }

  ## change tip labels:
  subtree$tip.label <- gsub("_", " ", (stringr::str_sub(subtree$tip.label, 17, -1)))

  species_data <- species_output |>
    dplyr::left_join(
      bacterial_taxonomy,
      by = join_by(genus),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      family = dplyr::if_else(genus %in% c("Nitrospira", "Spirochaeta"), NA_character_, family),
      order  = dplyr::if_else(genus %in% c("Spirochaeta"), NA_character_, order)
    ) |>
    dplyr::distinct()

  if (!is.null(genus)) species_to_include <- dplyr::filter(species_data, genus == .env$genus)
  if (!is.null(family)) species_to_include <- dplyr::filter(species_data, family == .env$family)
  if (!is.null(order)) species_to_include <- dplyr::filter(species_data, order == .env$order)
  if (!is.null(class)) species_to_include <- dplyr::filter(species_data, class == .env$class)

  species_to_include <- species_to_include |>
    dplyr::pull(species) |>
    unique()

  subtree <- castor::get_subtree_with_tips(subtree, only_tips = species_to_include)[[1]]

  data_tree <- tibble::as_tibble(subtree) |>
    dplyr::left_join(species_data, by = dplyr::join_by(label == species)) |>
    dplyr::filter(!is.na(node), !is.na(parent), !is.na(label)) |>
    tidyr::replace_na(list(
      family = "undefined",
      order  = "undefined",
      class  = "undefined",
      major_clade = "none"
    ))

  ## make the final tree data as a S4 object for plotting
  subtree_data <- tidytree::as.treedata(data_tree)

  max_dist_root <- max(get_all_distances_to_root(subtree))
  message("max_dist_root = ", max_dist_root)

  ## Plot tree
  p_tree <- ggtree::ggtree(
    subtree_data,
    layout = "rectangular",
    size = tree_line_size
  ) +
    ggtree::geom_tiplab(size = tip_text_size, fontface = "italic") +
    ggplot2::ggtitle(genus) +
    ggplot2::theme_classic(base_size = base_text_size) +
    ggplot2::theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.title = ggplot2::element_text(size = title_text_size, face = "bold"),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_text_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.title = ggplot2::element_text(size = legend_text_size),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  ## heatmap 1: resistance
  p1 <- ggtree::gheatmap(
    p_tree,
    species_data |> dplyr::select(species, resistance) |> tibble::column_to_rownames(var = "species"),
    offset = max_dist_root / 4,
    width = heatmap_width,
    color = NA,
    colnames = FALSE
  ) +
    ggplot2::scale_fill_manual(
      values = c("resistant" = "#ba181b", "susceptible" = "#d3d3d3"),
      labels = c("resistant", "susceptible", ""),
      na.value = hsv(0, 0, 0, 0),
      name = "Predicted\nresistance"
    )

  ## heatmap 2: evolvability
  p2 <- ggtree::gheatmap(
    p1 + ggnewscale::new_scale_fill(),
    species_data |> dplyr::select(species, evolvabilityI) |> tibble::column_to_rownames(var = "species"),
    offset = max_dist_root / 3,
    width = heatmap_width,
    color = NA,
    colnames = FALSE
  ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#0077b6", "#ffd700"),
      na.value = "grey20",
      name = "Predicted\nevolvability"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.title = ggplot2::element_text(size = legend_text_size)
    )

  fig_height <- max(6, length(subtree$tip.label) * tip_height)

  ggplot2::ggsave(
    filename = file_name,
    plot = p2,
    width = fig_width,
    height = fig_height,
    units = "in",
    limitsize = FALSE
  )

  return(p2)
}


#' plot subtree for a group of genera
#'
#' @param subtree a phylo object of the subset tree
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpsL
#'  sequences across bacterial species
#' @param bacterial_taxonomy a data frame providing info on bacterial taxonomic phylogeny
#' @param genera genera for which to plot a phylogenetic tree
#' @param families genera for which to plot a phylogenetic tree
#' @param orders genera for which to plot a phylogenetic tree
#' @param classes genera for which to plot a phylogenetic tree
#' @param file_name the path that the plot should be save in
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#'
plot_subtree_clades <- function(subtree, species_output, bacterial_taxonomy,
                                genera = NULL,
                                families = NULL,
                                orders = NULL,
                                classes = NULL,
                                file_path) {
  if (!is.null(genera)) {
    purrr::map(genera,
      \(x) plot_subtree_clade(
        subtree = subtree,
        species_output = species_output,
        bacterial_taxonomy = bacterial_taxonomy,
        genus = x,
        file_name = paste0(file_path, "/", x, ".pdf")
      ),
      .progress = TRUE
    )
  }
  if (!is.null(families)) {
    purrr::map(families,
      \(x) plot_subtree_clade(
        subtree = subtree,
        species_output = species_output,
        bacterial_taxonomy = bacterial_taxonomy,
        family = x,
        file_name = paste0(file_path, "/", x, ".pdf")
      ),
      .progress = TRUE
    )
  }
  if (!is.null(orders)) {
    purrr::map(orders,
      \(x) plot_subtree_clade(
        subtree = subtree,
        species_output = species_output,
        bacterial_taxonomy = bacterial_taxonomy,
        order = x,
        file_name = paste0(file_path, "/", x, ".pdf")
      ),
      .progress = TRUE
    )
  }
  if (!is.null(classes)) {
    purrr::map(classes,
      \(x) plot_subtree_clade(
        subtree = subtree,
        species_output = species_output,
        bacterial_taxonomy = bacterial_taxonomy,
        class = x,
        file_name = paste0(file_path, "/", x, ".pdf")
      ),
      .progress = TRUE
    )
  }
  return(invisible(NULL))
}

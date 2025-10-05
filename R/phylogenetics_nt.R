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

  
  return(list(permtest_resistance = permtest_resistance))
}

#' plot subtree nucleotide
#'
#' @param subtree a phylo object of the subset tree  
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpsL sequences across bacterial species
#' @param bacterial_taxonomy a data frame providing info on bacterial taxonomic phylogeny from NCBI
#' @param file_name the path that the plot should be save in
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#' @export
#'
#' @examples plot_subtree(subtree, filtered_output, meta_data, "./plots/myfilename.pdf")
plot_subtree_nt <- function(subtree, species_output, bacterial_taxonomy, file_name){
  
  # clades to be labeled in the tree:
  clades_to_label <- c(
    "Mycobacteriaceae",
    "Microbacteriaceae",
    # "Streptomycetaceae",
    # "Anaplasmataceae",
    "Lactobacillaceae",
    "Flavobacteriales", 
    # "Burkholderiaceae",
    "Burkholderiales",
    "Acidobacteriaceae", 
    # "Bacillaceae",
    # "Cyclobacteriaceae",
    "Cytophagales",
    "Streptosporangiales",
    "Clostridia",
    "Bacillales"
  )
  
  # change tip labels to make manipulations easier:  
  subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))
  
  species_data <- species_output |>
    left_join(bacterial_taxonomy, by = join_by(genus), relationship = "many-to-many") |> #join bacterial families to tree information by column species
    # #join with genus_variants to fix genus name
    # left_join(genus_variants |> select(genus_origin, family_var = family, order_var = order, class_var = class, phylum_var = phylum),
    #           by = c("genus" = "genus_origin"), relationship = "many-to-many") |> 
    # mutate(
    #   family = if_else(is.na(family), family_var, family),
    #   order = if_else(is.na(order), order_var, order),
    #   class = if_else(is.na(class), class_var, class),
    #   phylum = if_else(is.na(phylum), phylum_var, phylum)
    # ) |> 
    # select(-family_var, -order_var, -class_var, -phylum_var) |> 
    # There are multiple variants of these genera with different family and order assignments,
    # but the plot does not support duplicate species entries.    
    mutate(
      family = if_else(genus %in% c("Nitrospira", "Spirochaeta"), NA_character_, family),
      order  = if_else(genus %in% c("Spirochaeta"), NA_character_, order)
    ) |> 
    distinct() |>
    mutate(major_clade = NA) |>
    mutate(major_clade = ifelse(phylum %in% clades_to_label, phylum, major_clade)) |>
    mutate(major_clade = ifelse(class %in% clades_to_label, class, major_clade)) |>
    mutate(major_clade = ifelse(order %in% clades_to_label, order, major_clade)) |>
    mutate(major_clade = ifelse(family %in% clades_to_label, family, major_clade)) 
  # mutate(major_clade = ifelse(major_clade %in% c("Mollicutes", "Erysipelotrichales"),
  #                             "Mollicutes &\nErysipelotrichales", major_clade))
  
  data_tree <- as_tibble(subtree) |>
    left_join(species_data, by = join_by(label == species)) |>
    filter(!is.na(node), !is.na(parent), !is.na(label)) |>
    replace_na(list(family = 'undefined', 
                    order = 'undefined',
                    class = "undefined",
                    major_clade = "none"))
  
  ##identify common ancestor nodes for species in major bacterial orders
  clades <- data.frame(major_clade=unique(data_tree$major_clade),
                       common_ancestor=NA) |>
    filter(major_clade != "none")
  for (i in 1:length(clades$major_clade)) {
    clades$common_ancestor[i] <- getMRCA(subtree, data_tree$node[data_tree$major_clade == clades$major_clade[i]])
  }
  
  ##add clade nodes to data_tree
  data_tree <- full_join(data_tree, clades, by=join_by(major_clade))
  
  ##make the final tree data as a S4 object for plotting
  subtree_data <- tidytree::as.treedata(data_tree)
  
  ##plotting the tree
  p_tree <- ggtree(subtree_data, 
                   layout = "circular",
                   branch.length="none", 
                   color="#eee0cd")
  
  # add predicted resistance:  
  p1 <- gheatmap(p_tree, species_data |> select(species, resistance) |> column_to_rownames(var="species"), 
                 offset = 0, color=NULL, colnames=FALSE, width = 0.05) +
    scale_fill_manual(values=c("resistant" = "#ba181b", "susceptible" = "#d3d3d3"),
                      labels=c("resistant", "susceptible", ""),
                      na.value = hsv(0,0,0,0),
                      name = "Predicted\nresistance")
  
  # add major orders:
  #cols <- rep(brewer.pal(8, "Set1"), 100)[1:length(clades$major_clade)]
  cols <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Pastel1"))[1:length(clades$major_clade)]
  p2 <- p1 + geom_cladelab(node = clades$common_ancestor,
                           label = clades$major_clade,
                           align=TRUE,
                           geom='text',
                           #fill=cols,
                           fontsize = 3.5,
                           barcolour = cols, #"grey40",
                           offset.text = 2.5 ,
                           offset = 6,
                           barsize= 2, 
                           horizontal=FALSE,
                           hjust=0.5) +
    #angle = "auto") +
    theme(legend.position="bottom",
          legend.spacing.y = unit(0.1, "cm"))
  
  ggsave(filename = file_name, p2, width = 10, height = 10)
}

#' plot subtree for a specific clade (genus, family, order, class)
#'
#' @param subtree a phylo object of the subset tree  
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpsL
#'  sequences across bacterial species
#' @param bacterial_taxonomy a data frame providing info on bacterial taxonomic phylogeny from GTDB
#' @param genus_variants a data frame identifying genera in the filtered_output that correspond to multiple entries in the GTDB taxonomy data (e.g., "Actinomadura" represented as "Actinomadura_C" and "Actinomadura_D").
#' @param genus the name of the genus
#' @param family the name of the family
#' @param order the name of the order
#' @param class the name of the class
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#' @export
#'
#' @examples plot_subtree_nt(subtree, filtered_output, bacterial_taxonomy, "./plots/myfilename.pdf")
plot_subtree_clade_nt <- function(subtree, 
                               species_output, 
                               bacterial_taxonomy, 
                               genus = NULL, family = NULL, order = NULL, class = NULL,
                               file_name){
  
  if (sum(!is.null(c(genus, family, order, class))) != 1L)
    stop("One and only one of genus, family, order or class must be specified, the others need to be NULL.")
  
  ##change tip labels:  
  subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))
  
  species_data <- species_output |>
    left_join(bacterial_taxonomy, by = join_by(genus), relationship = "many-to-many") |> #join bacterial families to tree information by column species
    # There are multiple variants of these genera with different family and order assignments,
    # but the plot does not support duplicate species entries.    
    mutate(
      family = if_else(genus %in% c("Nitrospira", "Spirochaeta"), NA_character_, family),
      order  = if_else(genus %in% c("Spirochaeta"), NA_character_, order)
    ) |> 
    distinct() 
  
  if (!is.null(genus)) {
    species_to_include <- dplyr::filter(species_data, genus == .env$genus)
  }
  if (!is.null(family)) {
    species_to_include <- dplyr::filter(species_data, family == .env$family)
  }
  if (!is.null(order)) {
    species_to_include <- dplyr::filter(species_data, order == .env$order)
  }
  if (!is.null(class)) {
    species_to_include <- dplyr::filter(species_data, class == .env$class)
  }
  species_to_include <- species_to_include |>
    pull(species) |>
    unique()
  
  subtree <- castor::get_subtree_with_tips(subtree, only_tips = species_to_include)[[1]]
  
  data_tree <- as_tibble(subtree) |>
    left_join(species_data, by = join_by(label == species)) |>
    filter(!is.na(node), !is.na(parent), !is.na(label)) |>
    replace_na(list(family = 'undefined', 
                    order = 'undefined',
                    class = "undefined",
                    major_clade = "none"))
  
  ##make the final tree data as a S4 object for plotting
  subtree_data <- tidytree::as.treedata(data_tree)
  
  max_dist_root <- max(get_all_distances_to_root(subtree))
  print(max_dist_root)
  ##plotting the tree
  p_tree <- ggtree(subtree_data, 
                   layout = "rectangular") +
    geom_tiplab(size = 2) +
    ggtitle(genus)
  
  # add predicted resistance:  
  p1 <- gheatmap(p_tree, species_data |> select(species, resistance) |> column_to_rownames(var="species"), 
                 offset = max_dist_root/10, 
                 color=NULL, colnames=FALSE, 
                 width = 0.01) +
    scale_fill_manual(values=c("resistant" = "red", "susceptible" = "grey90"),
                      labels=c("resistant", "susceptible", ""),
                      na.value = hsv(0,0,0,0),
                      name = "Predicted\nresistance")
  
  # # add predicted evolvability:
  # p2 <- gheatmap(p1 + new_scale_fill(),
  #                species_data |> select(species, evolvabilityI) |> column_to_rownames(var="species"),
  #                offset = max_dist_root/9, color=NULL, colnames=FALSE, width = 0.01) +
  #   scale_fill_viridis_c(na.value = "grey20", name = "Predicted\nevolvability", option = "plasma") +
  #   scale_y_continuous(limits = c(0, round(1.005 * (nrow(data_tree) + 1) / 2))) +
  #   theme(legend.position="bottom")
  
  ggsave(filename = file_name, p1, 
         width = 15, 
         height = length(subtree$tip.label)/10,
         limitsize = FALSE)
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
plot_subtree_clades_nt <- function(subtree, species_output, bacterial_taxonomy,
                                genera = NULL, 
                                families = NULL, 
                                orders = NULL, 
                                classes = NULL, 
                                file_path){
  if (!is.null(genera)) {
    purrr::map(genera, 
               \(x) plot_subtree_clade_nt(subtree = subtree, 
                                       species_output = species_output, 
                                       bacterial_taxonomy = bacterial_taxonomy, 
                                       genus = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }
  if (!is.null(families)) {
    purrr::map(families, 
               \(x) plot_subtree_clade_nt(subtree = subtree, 
                                       species_output = species_output, 
                                       bacterial_taxonomy = bacterial_taxonomy, 
                                       family = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }
  if (!is.null(orders)) {
    purrr::map(orders, 
               \(x) plot_subtree_clade_nt(subtree = subtree, 
                                       species_output = species_output, 
                                       bacterial_taxonomy = bacterial_taxonomy, 
                                       order = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }  
  if (!is.null(classes)) {
    purrr::map(classes, 
               \(x) plot_subtree_clade_nt(subtree = subtree, 
                                       species_output = species_output, 
                                       bacterial_taxonomy = bacterial_taxonomy,
                                       class = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }
  return(invisible(NULL))
}

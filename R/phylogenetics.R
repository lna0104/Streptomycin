
#' build a tree based on GTDB bacterial life tree using species names and accessions
#'
#' @param output a data frame containing the predicted mutations for all species
#' @param original_tree a phylo object providing the bacterial tree of life (from GTDB)
#' @param meta_data a data frame providing the information of species which are included in bacterial tree of life (from GTDB)
#'
#' @return a list providing the subset tree and a data frame presenting all tip labels
#' @export
#'
#' @examples get_subtree(filtered_output, original_tree, meta_data)
get_subtree <- function(output, original_tree, meta_data){

  our_species <- output |>
    select(species, accession_numbers) |>
    distinct()
 
  #1. create a column "species_accessions" in metadata
  meta_data <- meta_data %>% 
    select(accession, ncbi_genbank_assembly_accession, ncbi_organism_name) %>%
    mutate(species_name = gsub("\\[|\\]", "", ncbi_organism_name)) %>% 
    mutate(species_name = gsub("\\'|\\'", "", species_name)) %>% 
    mutate(species_name = str_replace_all(species_name, "Candidatus", "")) %>%
    mutate(species_name = gsub("^ ", "", species_name)) %>% 
    mutate(species_name = sub("^([^ ]+[ ]+[^ ]+).*$", "\\1", #keeps only the first two words of a character element
                              species_name)) %>%
    filter(accession %in% original_tree$tip.label) %>%
    mutate(acc_ID = paste(ncbi_genbank_assembly_accession, 
                          species_name, sep = "_"))
  
  #2. change tip labels from "accessions to "species_accessions"
  original_tree$tip.label <- unname(setNames(meta_data$acc_ID, meta_data$accession)[original_tree$tip.label])
  original_tree$tip.label <- paste0(sub("A", "F", original_tree$tip.label))
  
  #3. exploring our species in metadata
  
  accession_not_tree <- setdiff(our_species$accession_numbers, substr(original_tree$tip.label, -1, 15))

  species_accession_not_tree <- our_species %>% 
    filter(accession_numbers %in% accession_not_tree) %>%
    pull(species)
  
  #4. subset original tree based on tip labels
  
  #subset the original tree based on presence of species names or accessions in tips:
  subset_tip_indices_species_name <- sapply(species_accession_not_tree, function(subset_label) {
    which(grepl(subset_label, original_tree$tip.label, fixed = TRUE))[1] #checks the presence of each species names in all tip labels, and returns index of first one
  })
  
  #subset based on accessions
  subset_tip_indices_accession <- sapply(our_species$accession_numbers, function(acc) {
    which(grepl(acc, original_tree$tip.label, fixed = TRUE))[1] #checks the presence of each species names in all tip labels, and returns index of first one
  })
  
  #subset for missing species
  # subset_tip_indices_missing_names <- sapply(missing_names$names_original_tree, function(subset_label) {
  #   which(grepl(subset_label, original_tree$tip.label, fixed = TRUE))[1] #checks the presence of each species names in all tip labels, and returns index of first one
  # }) 
  
  #subset for outliers to be removed
  # subset_tip_indices_outliers <- if (!is.null(outliers)) {
  #   which(purrr::reduce(
  #     lapply(outliers$outliers, function(subset_label) {
  #       grepl(subset_label, original_tree$tip.label, fixed = TRUE)
  #     }), 
  #     `|`  # Logical OR to combine all matches
  #   ))
  # } else {
  #   NULL
  # }
  
  #sum of all subsets to be kept
  subset_tip_indices <- unique(c(subset_tip_indices_species_name, 
                                 subset_tip_indices_accession))
  subset_tip_indices <- subset_tip_indices[!is.na(subset_tip_indices)]
  
  #6. get the sub_tree
  subset_tree <- get_subtree_with_tips(original_tree, only_tips = subset_tip_indices, 
                                      omit_tips= NULL, force_keep_root = TRUE)$subtree
  
  
  ### I did some extra coding to find actual name of those tip labels which only have "genus name sp."
  #make a data frame using subset_tree$tip.label
  tip_labels <- data.frame(old_tip_label = subset_tree$tip.label) %>%
    mutate(species_phylo = sub(".*_", "", old_tip_label)) %>%
    mutate(acc = substr(old_tip_label, 1, 15)) %>% #separate accession numbers as a column
    left_join(our_species, by = join_by(acc == accession_numbers)) %>%
    select(acc, species_phylo, species) %>% 
    distinct() %>%
    mutate(species = ifelse(is.na(species), species_phylo, species)) %>%
    mutate(new_tip_label = paste0(acc, "_", species))

  #replace new names in tip.labels of tree
  subset_tree$tip.label <- tip_labels$new_tip_label
  return(list(tree = subset_tree,
              tips = tip_labels))
}



#' Generate a phylogenetic tree with bootstrapping values
#'
#' @param seqs a XStringset object 
#' @return an object of class matrix which represents the binary matrix derived from the phylogenetic tree.
#' and at the same time returns the bootstrapping values for each node
#' @export
#'
#' @examples get_tree_with_boots(rpsL_target_sequences)
get_tree_with_boots <- function(seqs){
  #keep long enough sequences 
  seqs <- seqs[lengths(seqs) >= 2000]
  #perform multiple sequence alignment
  aligned_sequences <- muscle::muscle(seqs)
  # Convert alignment to a phyDat object
  phy_data <- as.phyDat(aligned_sequences)
  
  #compute all possible types of a given model for a given alignment
  estimated_models <- modelTest(phy_data, 
                                model="K80",
                                multicore = TRUE,
                                mc.cores = detectCores()) #count number of system cores
  #Generate or infer a tree using maximum likelihood (ML)
  computed_data <- pml_bb(estimated_models, 
                          rearrangement = "stochastic")
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


get_phylosignals <- function(subtree, 
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
  
  return(list(permtest_resistance = permtest_resistance,
              K_evolvabilityI = K_evolvabilityI,
              K_evolvabilityII = K_evolvabilityII,
              lambda_evolvabilityI = lambda_evolvabilityI,
              lambda_evolvabilityII = lambda_evolvabilityII))
}



#' plot subtree
#'
#' @param subtree a phylo object of the subset tree  
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpsL sequences across bacterial species
#' @param gtdb_taxonomy a data frame providing info on bacterial taxonomic phylogeny from GDTB
#' @param genus_variants a data frame identifying genera in the filtered_output that correspond to multiple entries in the GTDB taxonomy data (e.g., "Actinomadura" represented as "Actinomadura_C" and "Actinomadura_D").
#' @param file_name the path that the plot should be save in
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#' @export
#'
#' @examples plot_subtree(subtree, filtered_output, meta_data, "./plots/myfilename.pdf")
plot_subtree <- function(subtree, species_output, gtdb_taxonomy, genus_variants, file_name){
  
  # clades to be labeled in the tree:
  clades_to_label <- c(
    # "Sphingomonadaceae", #family
    "Sphingomonadales", #order
    "Devosiaceae", #family
    # "Pirellulales", #order
    # "Rhizobiales", #order
    # "Thermodesulfobacteriaceae", #family
    # "Actinopolymorphaceae", #family
    # "Micromonosporaceae", #family
    # "Mycobacteriaceae", #family
    # "Streptomycetaceae", #family
    # "Coriobacteriales", #order
    # "Actinomycetales", #order
    "Rickettsiales", #order
    "Micromonosporaceae", #family
    "Coriobacteriia", #class
    "Planctomycetia" #class
    # "Coriobacteriia",
  )
  
  # change tip labels to make manipulations easier:  
  subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))
  
  species_data <- species_output |>
    left_join(gtdb_taxonomy, by = join_by(genus), relationship = "many-to-many") |> #join bacterial families to tree information by column species
    #join with genus_variants to fix genus name
    left_join(genus_variants |> select(genus_origin, family_var = family, order_var = order, class_var = class, phylum_var = phylum),
      by = c("genus" = "genus_origin"), relationship = "many-to-many") |> 
    mutate(
      family = if_else(is.na(family), family_var, family),
      order = if_else(is.na(order), order_var, order),
      class = if_else(is.na(class), class_var, class),
      phylum = if_else(is.na(phylum), phylum_var, phylum)
    ) |> 
    select(-family_var, -order_var, -class_var, -phylum_var) |> 
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
                   color="grey85")
  
  # add predicted resistance:  
  p1 <- gheatmap(p_tree, species_data |> select(species, resistance) |> column_to_rownames(var="species"), 
                 offset = 0, color=NULL, colnames=FALSE, width = 0.05) +
    scale_fill_manual(values=c("resistant" = "red", "susceptible" = "grey90"),
                      labels=c("resistant", "susceptible", ""),
                      na.value = hsv(0,0,0,0),
                      name = "Predicted\nresistance")
  
  # add predicted evolvability:
  p2 <- gheatmap(p1 + new_scale_fill(),
                 species_data |> select(species, evolvabilityI) |> column_to_rownames(var="species"),
                 offset = 3.4, color=NULL, colnames=FALSE, width = 0.05) +
    scale_fill_viridis_c(na.value = "grey20", name = "Predicted\nevolvability", option = "plasma") +
    scale_y_continuous(limits = c(0, round(1.005 * (nrow(data_tree) + 1) / 2)))
  
  # add major orders:
  #cols <- rep(brewer.pal(8, "Set1"), 100)[1:length(clades$major_clade)]
  cols <- c(brewer.pal(9, "Set1"), brewer.pal(9, "Pastel1"))[1:length(clades$major_clade)]
  p3 <- p2 + geom_cladelab(node = clades$common_ancestor,
                           label = clades$major_clade,
                           align=TRUE,
                           geom='text',
                           #fill=cols,
                           fontsize = 3.5,
                           barcolour = cols, #"grey40",
                           offset.text = 2.5 ,
                           offset = 9,
                           barsize= 2, 
                           horizontal=FALSE,
                           hjust=0.5) +
                           #angle = "auto") +
    theme(legend.position="bottom")
  
  ggsave(filename = file_name, p3, width = 10, height = 10)
  
  # alternative labels using heatmap:
  # p4 <- gheatmap(p2 + new_scale_fill(),
  #                species_data |> select(species, major_clade) |> column_to_rownames(var="species"),
  #                offset = 9, color=NULL, colnames=FALSE, width = 0.05) +
  #   scale_fill_manual(values = cols, na.value = hsv(0,0,0,0), name = "Order")
  # ggsave(filename = "./plots/whole_genome_tree_test.pdf", p4, width = 15, height = 20)
}




#' plot subtree for a specific clade (genus, family, order, class)
#'
#' @param subtree a phylo object of the subset tree  
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpsL
#'  sequences across bacterial species
#' @param gtdb_taxonomy a data frame providing info on bacterial taxonomic phylogeny from GTDB
#' @param genus_variants a data frame identifying genera in the filtered_output that correspond to multiple entries in the GTDB taxonomy data (e.g., "Actinomadura" represented as "Actinomadura_C" and "Actinomadura_D").
#' @param genus the name of the genus
#' @param family the name of the family
#' @param order the name of the order
#' @param class the name of the class
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#' @export
#'
#' @examples plot_subtree(subtree, filtered_output, gtdb_taxonomy, "./plots/myfilename.pdf")
plot_subtree_clade <- function(subtree, 
                               species_output, 
                               gtdb_taxonomy, 
                               genus_variants,
                               genus = NULL, family = NULL, order = NULL, class = NULL,
                               file_name){
  
  if (sum(!is.null(c(genus, family, order, class))) != 1L)
    stop("One and only one of genus, family, order or class must be specified, the others need to be NULL.")
  
  ##change tip labels:  
  subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))
  
  species_data <- species_output |>
    left_join(gtdb_taxonomy, by = join_by(genus), relationship = "many-to-many") |> #join bacterial families to tree information by column species
    # Join with genus_variants to correct genus names split into variants 
    left_join(genus_variants |> select(genus_origin, family_var = family, order_var = order, class_var = class, phylum_var = phylum),
      by = c("genus" = "genus_origin"), relationship = "many-to-many") |> 
    mutate(
      family = if_else(is.na(family), family_var, family),
      order = if_else(is.na(order), order_var, order),
      class = if_else(is.na(class), class_var, class),
      phylum = if_else(is.na(phylum), phylum_var, phylum)
    ) |> 
    select(-family_var, -order_var, -class_var, -phylum_var) |> 
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
  
  # add predicted evolvability:
  p2 <- gheatmap(p1 + new_scale_fill(),
                 species_data |> select(species, evolvabilityI) |> column_to_rownames(var="species"),
                 offset = max_dist_root/9, color=NULL, colnames=FALSE, width = 0.01) +
    scale_fill_viridis_c(na.value = "grey20", name = "Predicted\nevolvability", option = "plasma") +
    scale_y_continuous(limits = c(0, round(1.005 * (nrow(data_tree) + 1) / 2))) +
    theme(legend.position="bottom")
  
  ggsave(filename = file_name, p2, 
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
#' @param gtdb_taxonomy a data frame providing info on bacterial taxonomic phylogeny from GTDB
#' @param genus_variants a data frame identifying genera in the filtered_output that correspond to multiple entries in the GTDB taxonomy data (e.g., "Actinomadura" represented as "Actinomadura_C" and "Actinomadura_D").
#' @param genera genera for which to plot a phylogenetic tree
#' @param families genera for which to plot a phylogenetic tree
#' @param orders genera for which to plot a phylogenetic tree
#' @param classes genera for which to plot a phylogenetic tree
#' @param file_name the path that the plot should be save in
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#'
plot_subtree_clades <- function(subtree, species_output, gtdb_taxonomy, genus_variants,
                                genera = NULL, 
                                families = NULL, 
                                orders = NULL, 
                                classes = NULL, 
                                file_path){
  if (!is.null(genera)) {
    purrr::map(genera, 
               \(x) plot_subtree_clade(subtree = subtree, 
                                       species_output = species_output, 
                                       gtdb_taxonomy = gtdb_taxonomy, 
                                       genus_variants = genus_variants,
                                       genus = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }
  if (!is.null(families)) {
    purrr::map(families, 
               \(x) plot_subtree_clade(subtree = subtree, 
                                       species_output = species_output, 
                                       gtdb_taxonomy = gtdb_taxonomy, 
                                       genus_variants = genus_variants,                                       
                                       family = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }
  if (!is.null(orders)) {
    purrr::map(orders, 
               \(x) plot_subtree_clade(subtree = subtree, 
                                       species_output = species_output, 
                                       gtdb_taxonomy = gtdb_taxonomy, 
                                       genus_variants = genus_variants,
                                       order = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }  
  if (!is.null(classes)) {
    purrr::map(classes, 
               \(x) plot_subtree_clade(subtree = subtree, 
                                       species_output = species_output, 
                                       gtdb_taxonomy = gtdb_taxonomy,
                                       genus_variants = genus_variants,
                                       class = x,
                                       file_name = paste0(file_path, "/", x, ".pdf")),
               .progress = TRUE)
  }
  return(invisible(NULL))
}


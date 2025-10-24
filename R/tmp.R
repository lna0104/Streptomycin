extract_gtdb_taxonomy <- function(filtered_output, meta_data){
  # Step 1: Extract unique species and accessions
  our_species <- filtered_output |>
    select(species, accession_numbers) |>
    distinct()
  
  meta_data_clean <- meta_data %>% 
    select(accession, ncbi_genbank_assembly_accession, ncbi_organism_name, gtdb_taxonomy) %>%
    mutate(ncbi_species_name = gsub("\\[|\\]", "", ncbi_organism_name)) %>% 
    mutate(ncbi_species_name = gsub("\\'|\\'", "", ncbi_species_name)) %>% 
    mutate(ncbi_species_name = str_replace_all(ncbi_species_name, "Candidatus", "")) %>%
    mutate(ncbi_species_name = gsub("^ ", "", ncbi_species_name)) %>% 
    mutate(ncbi_species_name = sub("^([^ ]+[ ]+[^ ]+).*$", "\\1", #keeps only the first two words of a character element
                                   ncbi_species_name)) %>%
    separate(gtdb_taxonomy,
             into = c("g_domain", "g_phylum", "g_class", "g_order", "g_family", "g_genus", "g_species"),
             sep = ";",
             remove = FALSE)  %>%
    mutate(across(
      c(g_domain, g_phylum, g_class, g_order, g_family, g_genus, g_species),
      ~ .x %>%
        str_replace("^[a-z]__","") %>%  # remove rank prefix
        na_if("") %>%                   # empty -> NA
        str_squish()
    )) %>%
    mutate(ncbi_genbank_assembly_accession = str_replace(ncbi_genbank_assembly_accession, "^GCA_", "GCF_")) %>%
    select(accession, ncbi_genbank_assembly_accession, ncbi_species_name, 
           g_domain, g_phylum, g_class, g_order, g_family, g_genus, g_species)
  
 
  
  # Step 3: Match species by accession
  species_by_accession <- our_species |>
    left_join(meta_data_clean, 
              by=c("accession_numbers" = "ncbi_genbank_assembly_accession")) 
  
  # Step 4: Handle unmatched species by name
  unmatched_species <- species_by_accession |>
    filter(is.na(accession))  |>
    pull(species)
  
  species_by_name <- our_species |>
    filter(species %in% unmatched_species) |>
    left_join(meta_data_clean, 
              by=c("species" = "ncbi_species_name")) 
  
  # Step 5: Combine matched species
  matched_species <- bind_rows(
    species_by_accession |>
      filter(!is.na(accession)),
    species_by_name)  
  
  final_species <- our_species |>
    left_join(matched_species, by = c("species", "accession_numbers")) |>
    distinct()
}

#' build a tree based on GTDB bacterial life tree using species names and accessions
#'
#' @param output a data frame containing the predicted mutations for all species
#' @param original_tree a phylo object providing the bacterial tree of life (from GTDB)
#' @param meta_data a data frame providing the information of species which are included in bacterial tree of life (from GTDB)
#'
#' @return 
#'   - tree: the subset phylogenetic tree
#' @export
#'
#' @examples extract_gtdb_subtree(filtered_output, original_tree, meta_data)
extract_gtdb_subtree <- function(original_tree, filtered_gtdb_taxonomy){
    
  
  # Step 6: Identify tip indices in the tree
  subset_tip_indices_accession <- sapply(filtered_gtdb_taxonomy$accession, function(acc) {
    which(grepl(acc, original_tree$tip.label, fixed = TRUE))[1] #checks the presence of each species names in all tip labels, and returns index of first one
  })
  
  #remove na 
  subset_tip_indices <- subset_tip_indices_accession[!is.na(subset_tip_indices_accession)]
  
  # Step 7: Extract subtree
  subset_tree <- get_subtree_with_tips(original_tree, only_tips = subset_tip_indices, 
                                       omit_tips= NULL, force_keep_root = TRUE)$subtree
  

  return(
    tree = subset_tree
  )
}



#' plot subtree
#'
#' @param subtree a phylo object of the subset tree  
#' @param filtered_output a data frame providing the results of mutation screening of high quality rpsL sequences across bacterial species
#' @param bacterial_taxonomy a data frame providing info on bacterial taxonomic phylogeny
#' @param file_name the path that the plot should be save in
#'
#' @return a plot presenting phylo_genetic relationship among bacterial species
#' @export
#'
#' @examples plot_subtree(subtree, filtered_output, meta_data, "./plots/myfilename.pdf")
plot_subtree_test <- function(subtree, species_output, filtered_gtdb_taxonomy, file_name){
  
  # clades to be labeled in the tree:
  clades_to_label <- c(
    "Sphingomonadales",
    "Devosiaceae",
    "Rickettsiales",
    "Coriobacteriia",
    "Planctomycetia",
    "Micromonosporaceae"
    # "Burkholderiales"
  )
  
  # change tip labels to make manipulations easier:  
  # subtree$tip.label <- gsub("_", " ", (str_sub(subtree$tip.label, 17, -1)))
  
  species_data <- species_output |>
    left_join(filtered_gtdb_taxonomy, by=join_by(species)) |>
    filter(!is.na(accession)) |>
    mutate(
      g_family = if_else(g_genus %in% c("Nitrospira", "Spirochaeta"), NA_character_, g_family),
      g_order  = if_else(g_genus %in% c("Spirochaeta"), NA_character_, g_order)
    ) |> 
    distinct() |>
    mutate(major_clade = NA) |>
    mutate(major_clade = ifelse(g_phylum %in% clades_to_label, g_phylum, major_clade)) |>
    mutate(major_clade = ifelse(g_class %in% clades_to_label, g_class, major_clade)) |>
    mutate(major_clade = ifelse(g_order %in% clades_to_label, g_order, major_clade)) |>
    mutate(major_clade = ifelse(g_family %in% clades_to_label, g_family, major_clade)) 
  # mutate(major_clade = ifelse(major_clade %in% c("Mollicutes", "Erysipelotrichales"),
  #                             "Mollicutes &\nErysipelotrichales", major_clade))
  

  data_tree <- as_tibble(subtree) |>
    left_join(species_data, by = join_by(label == accession)) |>
    filter(!is.na(node), !is.na(parent), !is.na(label)) |>
    replace_na(list(g_family = 'undefined', 
                    g_order = 'undefined',
                    g_class = "undefined",
                    major_clade = "none"))
  
  ##identify common ancestor nodes for species in major bacterial orders
  clades <- data.frame(major_clade=unique(data_tree$major_clade),
                       common_ancestor=NA) |>
    filter(major_clade != "none")
  for (i in 1:length(clades$major_clade)) {
    clades$common_ancestor[i] <- getMRCA(original_tree, data_tree$node[data_tree$major_clade == clades$major_clade[i]])
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
  p1 <- gheatmap(p_tree, species_data |> select(accession, resistance) |> column_to_rownames(var="accession"), 
                 offset = 0, color=NULL, colnames=FALSE, width = 0.05) +
    scale_fill_manual(values=c("resistant" = "#ba181b", "susceptible" = "#d3d3d3"),
                      labels=c("resistant", "susceptible", ""),
                      na.value = hsv(0,0,0,0),
                      name = "Predicted\nresistance")
  

  # add predicted evolvability:
  p2 <- p + gheatmap(
                 species_data |> select(accession, evolvabilityI) |> column_to_rownames(var="accession"),
                 offset = 3.4, color=NULL, colnames=FALSE, width = 0.05) +
    scale_fill_gradientn(colours = c("#0077b6", "#ffd700"), na.value = "grey20", name = "Predicted\nevolvability") +
    scale_y_continuous(limits = c(0, round(1.005 * (nrow(data_tree) + 1) / 2)))
  
  
  # add major orders:
  #cols <- rep(brewer.pal(8, "Set1"), 100)[1:length(clades$major_clade)]
  cols <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Pastel1"))[1:length(clades$major_clade)]
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

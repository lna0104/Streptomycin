#########################################################################
### Step 7: Phylogenetic distribution of resistance and evolvability  ###
#########################################################################

# input for this step:  filtered results for reliable sequences("./output/rpsL_filtered_output.csv")
#                       original bacterial phylogenetic tree of life ("./data/bac120.nwk") and 
#                       its metadata ("./data/bac120_metadata.tsv")
#                       NCBI bacterial taxonomic information ("./output/bacterial_taxonomy.csv")

# output for this step: subtree of original tree with tip_labels table ("./output/rpsL_subtree.RData")
#                       subtree of original tree (nwk file: "./output/rpsL_subtree.nwk")
#                       plot of subtree

# 1.load required files:
filtered_output <- read_csv("./output/rpsL_filtered_output.csv", show_col_types = FALSE)
original_tree <- read.tree("./data/bac120.nwk") #GTDB bacterial tree of life
# original_tree <- read.tree("./data/bac120.tree") #GTDB bacterial tree of life
bacterial_taxonomy <- read_csv("./data/rpsL_NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI
# meta_data <- read_tsv("./data/bac120_metadata.tsv", show_col_types = FALSE) #GTDB information on included species
gtdb_taxonomy <- read_csv("./data/gtdb_taxonomy.csv", show_col_types = FALSE)

# 2. get species-level summary of mutation screen data:
species_output <- get_species_output(filtered_output)

# 3.subset the tree based on species accessions and names:
subtree <- get_subtree(filtered_output, original_tree, meta_data)
write.tree(subtree$tree, file = "./output/rpsL_subtree.nwk") 

# 4. subtree visualization:
subtree <- read.tree("./output/rspL_subtree.nwk")
# big tree of all species:
plot_subtree(subtree, species_output, bacterial_taxonomy, file_name = "./plots/rpsL_phylogenies/whole_genome_tree.svg")
# smaller trees of individual clades:
plot_subtree_clades(subtree, species_output, bacterial_taxonomy, 
                    genera = c("Sphingomonas"),
                    families = c("Devosiaceae", "Mycobacteriaceae"),
                    orders = c("Pirellulales", "Sphingomonadales", "Rickettsiales"),
                    classes = c("Planctomycetia", "Alphaproteobacteria","Coriobacteriia"),
                    file_path = "./plots/rpsL_phylogenies/")

summarise_phylogenetics(subtree, species_output, sample_n = globsets$phylo_stats_sample_n, "./results/summary_rpsL_phylogenetics.txt")

#empty working environment to keep everything clean:
rm.all.but("globsets")


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
get_subtree <- function(output, original_tree, gtdb_taxonomy){
  
  our_species <- output |>
    select(species, accession_numbers) |>
    distinct()
  
  #1. create a column "species_accessions" in metadata
  meta_data_clean <- meta_data %>% 
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
  
  #2. change tip labels from "accessions" to "species_accessions"
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



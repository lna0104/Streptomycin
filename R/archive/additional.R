

##amino acids/codon lists

bases <- c("A", "C", "G", "T")
codons <- expand.grid(p1 = bases, p2 = bases, p3 = bases)
codons$codon <- paste(codons$p1, codons$p2, codons$p3, sep = "")
codons$aa <- NA

for(i in 1:nrow(codons)){
  DNAString(codons$codon[i])
  codons$aa[i] <- as.character(translate(DNAString(codons$codon[i])))
  
}
codons

##get all possible nucleotide compositions for all 20 aminoacids
some_codons <- codons$codon
poss_codons <- list()

for (i in 1:length(some_codons)){
  poss_codons[[i]] <- mutant_codons(some_codons[i])
} 

names(poss_codons) <- some_codons
poss_codons
poss_codons <- as.data.frame(poss_codons)
view(poss_codons)
poss_codons


which_mutations_in_Codon <- function(AA_pos, DNA_seq) {
  codon <- substr(DNA_seq, AA_pos*3 - 2, AA_pos*3)
  mut_codons <- mutant_codons(as.character(codon))
  return(mut_codons)
}

##pro alignment
# my_miss_match <- pairwiseAlignment(pattern = my_protein_ref, 
#                                  subject = my_protein_target)
# my_miss_match
# summary(my_miss_match)



# my_protein_target[511]
# my_protein_target[546]
#translate(substr(my_dna_target, 511*3 - 2, 511*3))
#translate(substr(my_dna_target, 511*3 - 2, 511*3))
#translate(my_dna_target)[511]


#bluh <- data.frame(pos = rep(c(1:1343), each = 3), nucleo = str_split_1(as.character(my_dna_target), ""), AA = rep(str_split_1(as.character(translate(my_dna_target)), ""), each = 3))

import_gff <- import(con = "./data/genomes/GCF_000006945.2_ASM694v2_genomic.gff",
                     format = "gff")

get_gene <- function(gbff_file, gene = "rpoB") {
  lines <- readLines(gbff_file)
  #We want to split the lines into a list of features.
  #Each feature begins with a line indented by exactly five(?!) spaces.
  features <- base::split(lines, cumsum(grepl("^     [^ ]", lines)))
  for(feature in features) {
    #We're only interested in CDSs.
    if(!grepl("^     CDS", feature[1])) {
      next
    }
    #Stitch all the lines together, stripping whitespace.
    text <- paste(trimws(feature), collapse = "")
    #Check if this is the right gene.
    if(!grepl(sprintf("\\gene=\"%s\"", gene), text)) {
      next
    }
    #Get the translation.
    translation <- gsub(".*/translation=\"([^\"]*)\".*", "\\1", text)
    return(translation)
  }
}


extract_species_names <- function(mutant){
    return(sub("[^_]*_([^_]*)_([^_]*)_.*", "\\1 \\2", mutant))
}


which_mutations_in_codon <- function(AA_pos, DNA_seq) {
  codon <- substr(DNA_seq, AA_pos*3 - 2, AA_pos*3)
  mut_codons <- mutant_codons(as.character(codon))
  return(mut_codons)
}



save_the_name <- function(rpoB_target_sequences){
  species_names <- c()
  species_names[i] <- summaries[[i]]$speciesname
  for (i in 1:length(summaries)){
    rpoB_target_sequences[i]
    paste0("<", species_names[i], sep = "")
  }
}  

which_AA <- function(AA_pos, DNA_seq) {
  AA <- translate(DNAString(substr(DNA_seq, AA_pos*3 - 2, AA_pos*3)))
  return(AA)
}


#plots
library(colorspace)
library(forcats)

num_mutation <- length(unique(mutations_list$AA_mut_name_Ecoli))
palette <- rwantshue::iwanthue()$hex(num_mutation)

mutations_list %>%
  subset(Species == c("Escherichia coli", "Mycobacterium tuberculosis")) %>%
  na.omit() %>%
  mutate(AA_pos_Ecoli = as.factor(AA_pos_Ecoli)) %>%
  ggplot(aes(AA_pos_Ecoli)) +
  geom_bar(aes(fill = AA_mut_name_Ecoli), width = 0.5) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.1)) + 
  labs(title ="The frequency of reported amino acid substituations in Escherichia coli and Mycobacterium tuberculosis", 
       x = "Amino acid positions", 
       y = "Number of unique entries") +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size = 8)) +
  theme(legend.key.size = unit(0.5, 'cm')) +
  guides(fill=guide_legend(title="Amino acid substituations")) +
  scale_fill_manual(values = palette)



#or
mutations_list %>%
  mutate(AA_mut_name_Ecoli = as.factor(AA_mut_name_Ecoli)) %>%
  na.omit() %>%
  ggplot() +
  geom_bar(aes(x = fct_infreq(AA_mut_name_Ecoli)), width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.1)) +
  labs(title ="The frequency of reported amino acid substituations across species",
       x = "Amino acid substituations",
       y = "Number of unique entries across species") +
  scale_y_continuous(breaks = seq(0, 30, by = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size = 0.2)) +
  theme(legend.key.size = unit(0.1, 'cm')) +
  guides(fill=guide_legend(title="Species")) +
  theme(legend.spacing = unit(0.1, "cm"))



#or
num_species <- length(unique(mutations_list$Species))
palette <- rwantshue::iwanthue()$hex(num_species)

mutations_list %>%
  mutate(AA_pos_Ecoli = as.factor(AA_pos_Ecoli)) %>%
  na.omit() %>%
  ggplot(aes(AA_pos_Ecoli)) +
  geom_bar(aes(fill= Species), width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.1)) + 
  labs(title ="The frequency of reported positions underwent amino acid substituation across species", 
       x = "Amino acid positions", 
       y = "Number of unique entries") +
  scale_y_continuous(breaks = seq(0, 120, by = 5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.2, 'cm')) +
  guides(fill=guide_legend(title="Species")) +
  theme(legend.spacing = unit(0.1, "cm")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = palette)


#or
mutations_list %>%
  mutate(AA_pos_Ecoli = as.factor(AA_pos_Ecoli)) %>%
  ggplot(aes(AA_pos_Ecoli)) +
  geom_bar(width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.1)) +
  labs(title ="The frequency of reported positions underwent amino acid substituation across species",
       x = "Amino acid positions",
       y = "Number of unique entries across species") +
  scale_y_continuous(breaks = seq(0, 110, by = 5)) +
  theme(plot.title = element_text(hjust = 0.5))


#for comparison: score of random sequence is around -10,000:
sample(c("A", "T", "G", "C"), length(rpoB_reference_Ecoli), replace = TRUE) |>
  paste0(collapse = "") |>
  DNAString() |>
  pairwiseAlignment(rpoB_reference_Ecoli) |>
  score()

sapply(str_split(filtered_final_output_df$species_names, " "), `[[`, 1)



#previous
frequency_mutation <-  muts %>%
  select(Species, AA_pos_Ecoli, AA_mutation, Origin) %>%
  filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) %>%
  mutate(mutation_name = paste(AA_pos_Ecoli, AA_mutation, sep = "_")) %>%
  distinct(Species, mutation_name) %>% #how many species have specific unique mutation 
  group_by(mutation_name) %>%
  filter(n() > 5) %>% #filter entries with frequency of < 5
  summarise(number_reports = n()) 


plot2 <- frequency_mutation %>%
  mutate(mutation_name = as.factor(mutation_name)) %>%
  ggplot() +
  geom_col(aes(x = reorder(mutation_name, -number_reports),  
               y = (number_reports)), 
           width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, 
                                   vjust=0.1, 
                                   hjust=0.95, 
                                   size=6), 
        axis.text.y = element_text(size=6), 
        axis.title=element_text(size=6, face="bold"),
        plot.title=element_text(hjust=0.5)) +
  labs(x = "Amino acid substitution",
       y = "Number of distinct entries across species") +
  scale_y_continuous(breaks = seq(0, 30, by=5), expand = c(0.01, 0)) +
  ggtitle('B')

#pie chart for mutation screen
mutation_screen_pai <- mutation_screen %>% 
  group_by(mutation_category) %>% 
  filter(n() > 0) %>% 
  summarise(number_reports = n()) 


mutation_screen_pai <- do.call(data.frame, mutation_screen_pai) %>% 
  mutate(perc = number_reports / sum(number_reports)) %>% 
  mutate(percentage = round(perc, digits = 4) * 100)  


plot7 <- mutation_screen_pai %>% 
  ggplot(aes(x="" , y=percentage, fill=mutation_category), position = "fill") +
  geom_col() +
  geom_text(aes(label = percentage), color = "white", size=3,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("grey20", "grey77", "red")) +
  labs(fill = "Mutation posibility:") +
  scale_y_continuous(labels = scales::percent)+
  coord_polar(theta = "y") +
  theme_classic() +
  theme_void() +
  theme(legend.position = "left")

plot7
#ggsave(filename = "./plots/mutation_screening_pie.pdf", plot7, width = 8, height = 6)


number_present_mutation2 <- plot5 + inset_element(
  plot4, 
  left = 0.1, 
  bottom = 0.1, 
  right = unit(1, 'npc') - unit(0.1, 'cm'), 
  top = unit(1, 'npc') - unit(0.1, 'cm')
)

plot5 <- plot4 + scale_y_continuous(limits = c(0, 300))


mutated_species <- filtered_final_output_df %>% 
  select(mutation_name, mutation_category, accession_numbers, genus_name) %>% 
  unique() %>% 
  filter(mutation_category == "present") %>% 
  group_by(mutation_name, genus_name) %>% 
  filter(n() > 2) %>% #filter out mutations with frequency of < 3 for each genus
  summarise(number_reports = n()) 


mutated_species <- mutated_species %>% 
  rbind(list('Others', 19)) 

mutated_species <- mutated_species %>% 
  mutate(perc = number_reports / sum(number_reports)) %>% 
  mutate(percentage = round(perc, digits = 4) * 100) 


plot7 <- mutated_species %>% 
  ggplot(aes(x="" , y=percentage, fill=mutation_name), position = "fill") +
  geom_col() +
  geom_text(aes(label = percentage), color = "black", size=2,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("olivedrab3", "mistyrose", "lightblue1", "darkslategray4", "khaki", "orange1")) +
  labs(fill = "Mutations:") +
  scale_y_continuous(labels = scales::percent)+
  coord_polar(theta = "y") +
  theme_classic() +
  theme_void() +
  theme(legend.position = "left",
        legend.text = element_text(size=6),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size=8)) 


plot7

highly_frequent_mutated_species <- plot6 + inset_element(plot7, left = 0.5, 
                                                         bottom = 0.3, 
                                                         right = unit(1, 'npc') - unit(0.1, 'cm'), 
                                                         top = unit(1, 'npc') - unit(0.1, 'cm'))
highly_frequent_mutated_species



# subset_labels <- str_remove(subset_labels, "CF_")
# subset_labels <- gsub("\\..*", "", subset_labels)

genomes <- list()
for (i in 1:length(genome_summeries)){
  species_names <- c()
  accession_numbers <- c()
  genomes$species_names[i] <- genome_summeries[[i]]$speciesname
  genomes$accession_numbers[i] <- genome_summeries[[i]]$assemblyaccession
  
}

###############################################


for (i in 1:length(rpoB_target_sequences_new)){
  print(is.null(rpoB_target_sequences_new[[i]]))
}


bluh <- as.data.frame(map(rpoB_target_sequences_new, is.null)) %>% 
  pivot_longer(cols = everything()) %>% 
  filter(value == TRUE)




duplicatedOnes <- rpoB_target_sequences[duplicated(names(rpoB_target_sequences))]
duplicatedOnes[sort(names(duplicatedOnes))]

taxonomy <- bacterial_taxonomy %>% 
  mutate(tax = str_replace_all(tax, "__", " ")) %>% 
  mutate(tax = str_remove(tax, gsub("^([^ ]+[ ]+[^ ]+[ ]+[^ ]+[ ]+[^ ]+[ ]+[^ ]+[ ]+).*$", "\\1", tax))) %>% 
  mutate(species = trimws(str_remove(tax, gsub("^([^ ]+[ ]+[^ ]+).*$", "\\1", tax)))) %>% 
  mutate(genus = sub(" .*", "", species)) %>% 
  mutate(family = trimws(gsub("^([^ ]+[ ]+).*$", "\\1", tax))) %>% 
  mutate(family = gsub("(.*);.*", "\\1", family)) %>% 
  select(species, genus, family) %>% 
  unique()


bacterial_taxonomy <- bacterial_taxonomy %>% 
  mutate(tax = str_replace_all(tax, "__", " ")) %>% 
  mutate(tax = trimws(str_remove(tax, gsub("^([^ ]+[ ]+[^ ]+[ ]+[^ ]+).*$", "\\1", tax)))) %>% 
  mutate(class = sub(" .*", "", tax)) %>% 
  mutate(class = gsub("(.*);.*", "\\1", class)) %>% 
  mutate(pre_order = trimws(str_remove(tax, sub(" .*", "", tax)))) %>% 
  mutate(order = sub(" .*", "", pre_order)) %>% 
  mutate(order = gsub("(.*);.*", "\\1", order)) %>% 
  mutate(family = sub(" .*", "", str_remove(pre_order, gsub("^([^ ]+[ ]+).*$", "\\1", pre_order)))) %>% 
  mutate(family = gsub("(.*);.*", "\\1", family)) %>% 
  select(class, order, family) %>% 
  unique()


bacterial_taxonomy <- full_join(bacterial_taxonomy, taxonomy, by=join_by(family))
write.csv(bacterial_taxonomy, file = "./output/bacterial_taxonomy.csv")


targets <- DNAStringSet(unlist(rpoB_target_sequences))
writeXStringSet(targets, filepath = "output/rpoB_target_sequences.fa")



#check what happend to species
final_output_checking <- final_output %>% select(species_names, accession_numbers) %>% 
  unique()
summaries_accessions <- c()
summaries_species <- c()
for (i in 1:length(genomes_info)){
  summaries_accessions[i] <- genomes_info[[i]]$assemblyaccession
  summaries_species [i] <- genomes_info[[i]]$speciesname
}

genomes_info <- data.frame(summaries_accessions, summaries_species)

summaries_accessions <- genomes_info$summaries_accessions
accession_numbers <- final_output_checking$accession_numbers

genome_database <- read_tsv("./data/ncbi_dataset.tsv")
missing <- genome_database %>% 
  filter(!genome_database_accesions %in% summaries_accessions) %>% 
  select(`Assembly Accession`, `Organism Name`)


get_alig_scores <- function(seqs, ref, region){
  ref_region <- ref[(region[1]*3 - 2):(region[2]*3)]
  ref_region_AA <- Biostrings::translate(ref_region)
  alig_scores <- data.frame(seq_name = names(seqs),
                            alig_score_region = NA)
  for(i in 1:nrow(alig_scores)) {
    # this probably won't work - need to find better function that ignores end gaps:
    alig <- pairwiseAlignment(pattern = Biostrings::translate(seqs[[i]]), subject = ref_region_AA, type = "global", )
  }
  return(alig_scores)
}


term <- 'GCF_000750005.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])'

no_seq <- sort(names(rpoB_target_sequences[lengths(rpoB_target_sequences) == 0]))

# #phylogenic analysis of sequences:
# high_quality_sequences <- unique(filtered_final_output_df$rpoB_target_sequence)
# 
# ##extracting the reliable sequences from main rpoB_target_sequences fasta file
# high_quality_sequences <- rpoB_target_sequences[names(rpoB_target_sequences) %in% long_enough_sequences]
# 
# ##save the sequences as a separate fasta file
# writeXStringSet(DNAStringSet(long_enough_sequences[!sapply(long_enough_sequences, is.null)]), #save rpoB sequences as a fasta file
#           filepath = "output/subset_sequences.fa") 


get_dist <- function(seqs, ref, start_pos, end_pos, maxdist){
  ref_region <- ref[(start_pos*3 - 2):(end_pos*3)]
  ref_region_AA <- Biostrings::translate(ref_region)
  RRDR_conservation_check <- data.frame(seq_name = extract_species_name(names(seqs)), 
                                        RRDR_conserved = NA, 
                                        RRDR_seq = NA, 
                                        distance = NA)
  for(i in 1:nrow(RRDR_conservation_check)) {
    con <- stringdist::extract(as.character(Biostrings::translate(seqs[[i]], if.fuzzy.codon = "X")), 
                               pattern = as.character(ref_region_AA), 
                               maxDist = maxdist, 
                               method = "lv")
    RRDR_conservation_check$RRDR_seq[i] <- ifelse(is.na(con[1, 1]), "None", con[1, 1])
    RRDR_conservation_check$RRDR_conserved[i] <- !is.na(con[1, 1])
    RRDR_conservation_check$distance[i] <- ifelse(is.na(con[1, 1]), "more than 50", 
                                                  stringdist(as.character(ref_region_AA), 
                                                             con[1, 1],
                                                             method = "lv"))}}

#another possible plots
# view_tree <- ggtree(subtree_data, aes(color=resistancy, size = resistancy, angle = angle), layout = "circular", branch.length="none")+
#   geom_tiplab(aes(label = ifelse(major_species == "TRUE", label, "")), size = 1, align = T) #show tip labels only for major species
# view_tree

# view_tree <- ggtree(subtree_data, color="grey85", layout = "circular", branch.length="none")+
#   geom_tippoint(aes(fill=resistancy, color=resistancy), size=1, shape=21) + #show tip labels only for major species
#   annotate("point", color="black")
# view_tree



# view_tree <- ggtree(subtree_data, layout = "circular", branch.length="none", color="grey85") +
#   geom_point2(aes(subset=(resistancy=="resistant")), shape=21, size=0.01, fill='red', color='red')+
#   scale_size_manual(values = c(0.1 ,0.1)) +
#   geom_cladelab(node=clades$common_ancestor[4], label="Enterobacteriaceae", align=TRUE, 
#                 geom='label', fill='lightskyblue', fontsize = 8, barcolour = 'lightskyblue', offset.text = 5 , offset = -1) +
#   geom_hilight(node=clades$common_ancestor[4], fill='lightskyblue', type="rect")+
#   geom_cladelab(node=clades$common_ancestor[2], label="Mycobacteriaceae", align=TRUE, 
#                 geom='label', fill='yellow', fontsize = 8, barcolour = 'yellow', offset.text = 5 , offset = -1) +
#   geom_hilight(node=clades$common_ancestor[2], fill='yellow', type="rect") +
#   geom_cladelab(node=clades$common_ancestor[3], label="Streptomycetaceae", align=TRUE, 
#                 geom='label', fill='mistyrose', fontsize = 8, barcolour = 'mistyrose', offset.text = 5 , offset = -1) +
#   geom_hilight(node=clades$common_ancestor[3], fill='mistyrose', type="rect") +
#   geom_cladelab(node=clades$common_ancestor[6], label="Pseudomonadaceae", align=TRUE, 
#                 geom='label', fill='darkgoldenrod2', fontsize = 8, barcolour = 'darkgoldenrod2', offset.text = 5 , offset = -1) +
#   geom_hilight(node=clades$common_ancestor[6], fill='darkgoldenrod2', type="rect") +
#   geom_cladelab(node=clades$common_ancestor[5], label="Moraxellaceae", align=TRUE, 
#                 geom='label', fill='darkolivegreen3', fontsize = 8, barcolour = 'darkolivegreen3', offset.text = 5 , offset = -1) +
#   geom_hilight(node=clades$common_ancestor[5], fill='darkolivegreen3', type="rect")
# view_tree


bacterial_taxonomy <- bacterial_taxonomy %>% 
  select(Lineage, `Scientific name`) %>% 
  mutate(family = str_extract(Lineage, "[A-Za-z]+ceae")) %>% 
  mutate(order = str_extract(Lineage, "[A-Za-z]+ales")) %>% 
  mutate(genus = sub(" .*", "", `Scientific name`)) %>% 
  mutate(genus = gsub("\\[|\\]", "", genus)) %>%
  mutate(genus = str_replace_all(genus, "'", "")) %>% 
  mutate(genus = sub("/.*", "", genus)) %>% 
  mutate(genus = gsub("^[a-z].*\\s?", "none", genus)) %>% 
  mutate(genus = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "none", genus)) %>% 
  mutate(family = ifelse(str_detect(genus, "[A-Za-z]+ceae"), genus, family)) %>% 
  mutate(order = ifelse(str_detect(genus, "[A-Za-z]+ales"), genus, order)) %>% 
  mutate(genus = str_replace_all(genus, c("ales" = "", 
                                          "ceae" = "", 
                                          "[A-Za-z]+-like" = "none",
                                          "[A-Za-z]+(Iran)" = "none"))) %>%
  filter(!genus %in% terms_to_remove) %>% 
  filter(genus != "none") %>% 
  filter(!is.na(family), !is.na(order)) %>% 
  filter(!(genus == "Corynebacterium" & family == "Erwiniaceae")) %>% 
  filter(!(genus == "Enterobacter" & family == "Erwiniaceae")) %>% 
  filter(!(genus == "Curtobacterium" & family == "Erwiniaceae")) %>% 
  select(order, family, genus) %>% 
  unique()


bacterial_taxonomy <- bacterial_taxonomy %>%
  select(family, genus) %>% 
  unique()


element_counts <- table(bacterial_taxonomy$genus)
redundant_elements <- names(element_counts[element_counts > 1])

bacterial_taxonomy <- bacterial_taxonomy %>%
  select(family, genus) %>%
  mutate(family = case_when(genus %in% redundant_elements ~ "undefined", TRUE ~ family)) %>% 
  distinct()


bacterial_taxonomy <- bacterial_taxonomy %>% 
  mutate(tax = str_replace_all(tax, "__", " ")) %>%
  mutate(tax = str_replace_all(tax, "-", "")) %>%
  mutate(tax = str_replace_all(tax, "_", "")) %>%
  mutate(family = sub(".*;f", "", tax)) %>% 
  mutate(family = sub(";.*", "", family)) %>% 
  mutate(family = str_replace_all(family, " ", " ")) %>%
  mutate(family = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "undefined", family)) %>% 
  mutate(family = sub("[A-Z]$", "", family)) %>% 
  mutate(order = sub(".*;o", "", tax)) %>% 
  mutate(order = sub(";.*", "", order)) %>% 
  mutate(order = str_replace_all(order, " ", " ")) %>% 
  mutate(order = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "undefined", order)) %>% 
  mutate(order = sub("[A-Z]$", "", order)) %>% 
  mutate(genus = sub(".*;g", "", tax)) %>% 
  filter(!is.na(family), !is.na(order)) %>% 
  mutate(genus = sub("^\\s+", "", genus)) %>% 
  mutate(genus_1 = sub(" .*", "", genus)) %>% 
  mutate(genus_1 = str_replace_all(genus_1, ";s", "")) %>%
  mutate(genus_1 = str_replace_all(genus_1, "_", "")) %>% 
  mutate(genus_1 = str_replace_all(genus_1, "-", "")) %>% 
  mutate(genus_1 = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "undefined", genus_1)) %>% 
  mutate(genus_1 = sub("[A-Z]$", "", genus_1)) %>% 
  select(genus_1, family, order) %>% 
  unique() %>% 
  mutate(genus = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "undefined", genus_1)) %>% 
  mutate(family = gsub("^[0-9]+$", "undefined", family)) %>% 
  mutate(order = gsub("^[0-9]+$", "undefined", order)) %>% 
  mutate(genus = gsub("^[0-9]+$", "undefined", genus)) %>% 
  select(genus, family, order) %>% 
  mutate(genus = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "undefined", genus)) %>% 
  mutate(order = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "undefined", order)) %>% 
  mutate(family = gsub("[a-zA-Z]+\\d+|\\d+[a-zA-Z]+", "undefined", family)) %>% 
  mutate(genus = gsub("\\d+", "undefined", genus)) %>% 
  mutate(genus = str_replace_all(genus, "undefinedundefined", "undefined")) %>% 
  mutate(order = str_replace_all(order, "undefinedundefined", "undefined")) %>%
  mutate(family = str_replace_all(family, "undefinedundefined", "undefined")) %>%
  filter(!grepl("undefined", genus)) %>% 
  filter(!grepl("undefined", order)) %>% 
  filter(!grepl("undefined", family)) %>% 
  mutate(genus = sub("^\\s+", "", genus)) %>% 
  mutate(family = sub("^\\s+", "", family)) %>% 
  mutate(order = sub("^\\s+", "", order)) 




new <- data.frame(genus = "Enterobacteriaceae",
                  family = "Enterobacteriaceae",
                  order = "Enterobacterales")

refined_taxonomy <- bacterial_taxonomy %>%
  mutate(genus = sub("[A-Z]$", "", genus)) %>% 
  mutate(order = case_when(genus == "Synechococcus" & order == "Neosynechococcales" ~ "Synechococcales",
                           TRUE ~ order)) %>% 
  mutate(family = case_when(genus == "Synechococcus" & family == "Thermosynechococcaceae" ~ "Synechococcaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Spirochaeta" & family == "Salinispiraceae" ~ "Spirochaetaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Thioalkalivibrio" & family == "Thioalkalivibrionaceae" ~ "Ectothiorhodospiraceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Pseudomonas" & family == "Xanthomonadaceae" ~ "Pseudomonadaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Pseudomonas" & order == "Xanthomonadales" ~ "Pseudomonadales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Rhodospirillum" & family == "Azospirillaceae" ~ "Rhodospirillaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Rhodospirillum" & order == "Azospirillales" ~ "Rhodospirillales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Ruminiclostridium" & family == "Ruminococcaceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Ruminiclostridium" & family == "Acetivibrionaceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Ruminiclostridium" & family == "Ethanoligenenaceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Ruminiclostridium" & order == "Acetivibrionales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(order = case_when(genus == "Ruminiclostridium" & order == "Oscillospirales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Ruminococcus" & family == "Ruminococcaceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Ruminococcus" & family == "Acutalibacteraceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Ruminococcus" & family == "Lachnospiraceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Ruminococcus" & order == "Oscillospirales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(order = case_when(genus == "Ruminococcus" & order == "Lachnospirales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Phormidesmis" & family == "Phormidesmiaceae" ~ "Leptolyngbyaceae",
                            TRUE ~ family)) %>% 
  mutate(order = case_when(genus == "Phormidesmis" & order == "Phormidesmiales" ~ "Leptolyngbyales",
                           TRUE ~ order)) %>% 
  mutate(family = case_when(genus == "Phormidium" & family == "Geitlerinemaceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Phormidium" & family == "Phormidiaceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Phormidium" & order == "Cyanobacteriales" ~ "Oscillatoriales",
                           TRUE ~ order)) %>%
  mutate(order = case_when(genus == "Phormidium" & order == "Burkholderiales" ~ "Polyangiales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Phormidium" & family == "Burkholderiaceae" ~ "Polyangiaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Prosthecomicrobium" & order == "Rhizobiales" ~ "Hyphomicrobiales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Prosthecomicrobium" & family == "Ancalomicrobiaceae" ~ "Kaistiaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Eubacterium" & order == "Lachnospirales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Eubacterium" & family == "Lachnospiraceae" ~ "Eubacteriaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Eubacterium" & order == "Peptostreptococcales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Eubacterium" & family == "Anaerovoracaceae" ~ "Eubacteriaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Eubacterium" & order == "Oscillospirales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Eubacterium" & family == "Acutalibacteraceae" ~ "Eubacteriaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Leptolyngbya" & order == "Elainellales" ~ "Leptolyngbyales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Leptolyngbya" & family == "Elainellaceae" ~ "Leptolyngbyaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Mycoplasma" & family == "Metamycoplasmataceae" ~ "Mycoplasmataceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Mycoplasma" & family == "Mycoplasmoidaceae" ~ "Mycoplasmataceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Paenibacillus" & family == "Reconcilibacillaceae" ~ "Paenibacillaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Bacteroides" & order == "Lachnospirales" ~ "Bacteroidales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Bacteroides" & family == "Lachnospiraceae" ~ "Bacteroidaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Clostridium" & order == "Clostridiales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Erysipelotrichaceae" ~ "Bacteroidaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Clostridium" & order == "Peptostreptococcales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(order = case_when(genus == "Clostridium" & order == "Tissierellales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Natronincolaceae" ~ "Clostridiaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Tepidimicrobiaceae" ~ "Clostridiaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Lachnospiraceae" ~ "Clostridiaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Caloramatoraceae" ~ "Clostridiaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Bacteroidaceae" ~ "Clostridiaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Acutalibacteraceae" ~ "Clostridiaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Clostridium" & family == "Peptostreptococcaceae" ~ "Clostridiaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Clostridium" & order == "Lachnospirales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Desulfotomaculum" & family == "Desulfovirgulaceae" ~ "Desulfotomaculaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Bacillus" & family == "Caldibacillaceae" ~ "Bacillaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Bacillus" & family == "Caldalkalibacillaceae" ~ "Bacillaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Bacillus" & family == "Salisediminibacteriaceae" ~ "Bacillaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Bacillus" & family == "Anaerobacillaceae" ~ "Bacillaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Bacillus" & family == "Marinococcaceae" ~ "Bacillaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Bacillus" & family == "Domibacillaceae" ~ "Bacillaceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Bacillus" & family == "Desulfovirgulaceae" ~ "Bacillaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Acetivibrio" & order == "Lachnospirales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Acetivibrio" & family == "Acetivibrionaceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Acetivibrio" & order == "Acetivibrionales" ~ "Eubacteriales",
                           TRUE ~ order)) %>%
  mutate(family = case_when(genus == "Acetivibrio" & family == "Lachnospiraceae" ~ "Oscillospiraceae",
                            TRUE ~ family)) %>%
  mutate(family = case_when(genus == "Polyangium" & family == "Burkholderiaceae" ~ "Polyangiaceae",
                            TRUE ~ family)) %>%
  mutate(order = case_when(genus == "Polyangium" & order == "Burkholderiales" ~ "Polyangiales",
                           TRUE ~ order)) %>% 
  unique()
refined_taxonomy <- rbind(new, refined_taxonomy)

write.csv(refined_taxonomy, "./data/bacterial_taxonomy_recent.csv")


data_tree <- as_tibble(sub_tree) %>% 
  mutate(genus = sub(" .*", "", label)) %>% #add a column for genus
  full_join(bacterial_taxonomy, relationship = "many-to-many") %>% #join bacterial families to tree information by column genus
  unique() %>% 
  filter(!is.na(node), !is.na(parent), !is.na(label)) %>% 
  replace_na(list(family = 'undefined', 
                  order = 'undefined')) %>%
  mutate(major_family = ifelse(family == "Mycobacteriaceae", "Mycobacteriaceae", #we need this column for highlighting some major bacterial families
                               ifelse(family == "Enterobacteriaceae", "Enterobacteriaceae", 
                                      ifelse(family == "Treponemataceae", "Treponemataceae",
                                             ifelse(family == "Pseudomonadaceae", "Pseudomonadaceae",
                                                    ifelse(family == "Streptomycetaceae", "Streptomycetaceae",
                                                           ifelse(family == "Streptococcaceae", "Streptococcaceae",
                                                                  ifelse(family == "Bifidobacteriaceae", "Bifidobacteriaceae",
                                                                         ifelse(family == "Leptospiraceae", "Leptospiraceae",
                                                                                ifelse(family == "Metamycoplasmataceae", "Metamycoplasmataceae",
                                                                                       ifelse(family == "Mycoplasmataceae", "Mycoplasmataceae",
                                                                                              ifelse(family == "Spiroplasmataceae", "Spiroplasmataceae",
                                                                                                     ifelse(family == "Microbacteriaceae", "Microbacteriaceae",
                                                                                                            ifelse(family == "Thermomonosporaceae", "Thermomonosporaceae",
                                                                                                                   ifelse(family == "Lactobacillaceae", "Lactobacillaceae",
                                                                                                                          ifelse(family == "Borreliaceae", "Borreliaceae",
                                                                                                                                 ifelse(family == "Moraxellaceae", "Moraxellaceae", 
                                                                                                                                        ifelse(family == "Corynebacteriaceae", "Corynebacteriaceae", 
                                                                                                                                               ifelse(family == "Streptosporangiaceae", "Streptosporangiaceae", 
                                                                                                                                                      ifelse(family == "Actinomycetaceae", "Actinomycetaceae",
                                                                                                                                                             ifelse(family == "Brachyspiraceae", "Brachyspiraceae",
                                                                                                                                                                    ifelse(family == "Erysipelotrichaceae", "Erysipelotrichaceae",
                                                                                                                                                                           ifelse(family == "Selenomonadaceae", "Selenomonadaceae",
                                                                                                                                                                                  ifelse(family == "Nocardiaceae", "Nocardiaceae",
                                                                                                                                                                                         "none")))))))))))))))))))))))) %>% 
  mutate(resistancy = ifelse(label %in% resistant_species, #we need this column to color the tree based on
                             "resistant", "sensitive")) %>% 
  mutate(major_species = ifelse(label %in% major_species, #we need this column to only show labels for highly frequent species in muts
                                "TRUE", "FALSE"))


##identify common ancestor nodes for species in major bacterial families
clades <- data.frame(major_family=unique(data_tree$major_family), 
                     common_ancestor=NA)

for (i in 1:length(clades$major_family)) {
  clades$common_ancestor[i] <- getMRCA(sub_tree, data_tree$node[data_tree$major_family == clades$major_family[i]]
  )
}



#Proportion (Fraction) Normalization
proportion_normalize <- function(x) {
  return(x / sum(x))
}

#Min-Max Normalization
min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

#z_score normalizization
# Calculate mean and standard deviation for normalization
mean_n_pos <- mean(num_poss$n_pos)
sd_n_pos <- sd(num_poss$n_pos)

# Add normalized column
num_poss <- num_poss %>%
  mutate(norm = (n_pos - mean_n_pos) / sd_n_pos)

# Shift normalized values to be positive
min_norm <- min(num_poss$norm)
num_poss <- num_poss %>%
  mutate(norm = norm + abs(min_norm))
##############################################################################
### Step 0: Loading all libraries and scripts and define global settings ###
##############################################################################

library(ALJEbinf)
library(tidyverse)
library(openssl)
library(stringr)
library(parallel)
library(future.apply)
library(varhandle)

library(Biostrings)
library(pwalign)
library(MSA2dist)

library(castor)
library(ape)
library(phytools)
library(phangorn)
library(ggtree)
library(tidytree)
library(treeio)
library(bio3d)

library(ggnewscale)
library(colorspace)
library(scales)
library(patchwork)
library(ggpubr)
library(plotrix)
library(RColorBrewer)
library(GGally)
library(ggh4x)
library(pander)
library(quarto)
library(NGLVieweR)
library(htmlwidgets)

source("R/util.R")
source("R/bioinformatics.R")
source("R/analyses.R")
source("R/plotting.R")
source("R/reports.R")
source("R/phylogenetics.R")
source("R/codon_networks.R")
#source("R/Yang_data.R")
source("R/structure.R")

# Global settings and parameters:
globsets <- list(
  min_n_studies = 3, # minimum number of studies that a mutation needs to be reported in for inclusion
  min_n_species = 3, # minimum number of species that a mutation needs to be reported in for inclusion
  min_seq_length = 300, # minimum length of included gene target sequences
  min_alig_score = -Inf, # minimum alignment score (with E. coli) of included gene target sequences
  max_core_dist = 35, # maximum Levenshtein distance between E. coli core gene region to corresponding target region
  phylo_stats_sample_n = 5000, # number of species to sample for phylogenetics statistics
  # Yang_min_log2_fc = 0, # minimum log2(foldchange) in frequencies to include a mutation as "selected for"
  random_seed = 22)
options(nwarnings = 10000)

set.seed(globsets$random_seed)


########################################################################
### Step 1: Processing table of known STR resistance mutations       ###
########################################################################

# input for this step: table of rpsL reference sequences ("./data/rpsL_references.csv")
#                      fasta file with rpsL reference sequences (""./data/rpsL_references.fasta"")
#                      table of reported STR mutations ("./data/reported_mutations.csv")
# output for this step: coordinates ("./output/coordinates.RData"), 
#                       table of reported mutations ("./output/muts.csv")
#                       a plot showing distributions of reported mutations ("./plots/reported_mutations.pdf)
#                       a summary of the reported mutations ("./results/summary_reported_mutations.txt")

# 1 get rpsL ref sequences coordination of all reference sequences against the one in Escherichia coli

#load reference sequences and their information 
refs <- read_csv("./data/rpsL_references.csv", show_col_types = FALSE)
seqs <- readDNAStringSet("./data/rpsL_references.fasta")
mutation_list_reports <- read_csv("./data/reported_mutations.csv", show_col_types = FALSE)|>
  filter(Gene=="rpsL")

# Check whether the above files have been changed and hence the coordinates need to be updated 
if(file.exists("./data/fastahash.Rds") && as.character(openssl::sha1(file("./data/rpsL_references.fasta"))) == readRDS("./data/fastahash.Rds")){
  print("Sequences file has not changed, loading original coordinates")
  load(file = "./output/coordinates.RData")
} else {
  old_fastahash <- as.character(openssl::sha1(file("./data/rpsL_references.fasta")))
  saveRDS(old_fastahash, 
       file = "./data/fastahash.Rds")
  print("Sequences file has changed, regenerating coordinates")
  #get coordinates
  coordinates <- ALJEbinf::getAllCoordinates(seqs, "rpsL_Escherichia_coli_MG1655")
  save(coordinates, file = "./output/coordinates.RData")
}

#2 load and complete table of reported mutations:
muts <- mutation_list_reports |>
  fillMutationsTable(refs, seqs, coordinates) |>
  filter(!is.na(AA_mut_name_Ecoli))  # filter out "bad" entries that couldn't be mapped to E. coli

write_csv(muts, "./output/muts.csv")

#3 summary and plot of reported mutations
plot_reported_mutations(muts, file_name = "./plots/reported_mutations.pdf", n_frequency = 3) # returns frequent reported mutations, positions and species
summarise_reported_mutations(muts, file_name = "./results/summary_reported_mutations.txt") # returns a text message summarizing previous reports

#empty working environment to keep everything clean
rm.all.but("globsets")

##############################################################################
### Step 2: Retrieving target sequences from all bacterial reference genomes ###
##############################################################################

# input for this step: none, apart from search term and database specified below
# output for this step: summaries of downloaded genomes ("./output/summaries.rds"), 
#                       genomic files ("./output/genomes/*.fna")
#                       extracted target sequences (example rpsL) ("./output/rpsL_target_sequences.fa")
#                       taxonomy database for all downloaded genera ("./output/NCBI_taxonomy.csv")

# define the desired database
db <- "assembly" 
# define search term for representative genomes: Bacterial genomes at all assembly levels with annotation
# term_rep <- '("Escherichia coli"[Organism] OR "Mycobacterium"[Organism] OR "Mycoplasmatota"[Organism]) AND ("latest refseq"[filter] AND "representative genome"[filter] AND "refseq has annotation"[Properties])'#representative

term_rep <- '("Bacteria"[Organism] OR bacteria[All Fields]) AND ("latest refseq"[filter] AND "representative genome"[filter] AND "refseq has annotation"[Properties])'#representative
#search the entire database using the defined term
summaries_rep <- get_summaries(db, term_rep)

# define search term for reference genomes: Bacterial genomes at all assembly levels with annotation 
# term_ref <- '("Escherichia coli"[Organism] OR "Mycobacterium"[Organism] OR "Mycoplasmatota"[Organism] AND ("latest refseq"[filter] AND "reference genome"[filter] AND "refseq has annotation"[Properties])'

term_ref <- '("Bacteria"[Organism] OR bacteria[All Fields]) AND ("latest refseq"[filter] AND "reference genome"[filter] AND "refseq has annotation"[Properties])'
# search the entire database using the defined term
summaries_ref <- get_summaries(db, term_ref)

# combine summaries and save them
summaries <- c(summaries_ref, summaries_rep)
saveRDS(summaries, file = "./output/summaries.rds")

# download genomes
download_files(summaries, dir = "output/genomes")

# extract target sequences and save them as rds and fasta
rpsL_target_sequences <- get_target_sequences(summaries, 
                                              dir = "output/genomes", 
                                              target_gene = 'rpsL', 
                                              target_protein = "30S ribosomal subunit protein S12|30S ribosomal protein S12")
writeXStringSet(DNAStringSet(rpsL_target_sequences), filepath = "output/rpsL_target_sequences.fa")

# download and extract taxonomy information for downloaded genomes
download_taxonomy(summaries, output_file = "./data/NCBI_taxonomy.csv")

#empty working environment to keep everything clean
rm.all.but("globsets")

#########################################################################
### Step 3: Checking all target sequences for reported mutations      ###
#########################################################################

# input for this step:  processed table of reported mutations ("./output/muts.csv")
#                       reported frequencies of rpsL mutations under STR 
#                       E. coli rspL reference sequence (from "./data/rspL_references.fasta")
#                       extracted rspL target sequences ("./output/rspL_target_sequences.fa")
# output for this step: Table of screening results, including existing and evolvable mutations 
#                           for all target sequences ("./output/raw_output.csv")

# 1.load required data:
muts <- read.csv("./output/muts.csv")
rpsL_target_sequences <- readDNAStringSet("./output/rpsL_target_sequences.fa")
rpsL_reference_Ecoli <- readDNAStringSet("./data/rpsL_references.fasta")[["rpsL_Escherichia_coli_MG1655"]]

# 2.make a list of reliable mutations to be screened:
mutation_list_reports <- filter_mutations(muts,
                                          min_n_species = globsets$min_n_species, 
                                          min_n_studies = globsets$min_n_studies)
mutation_list <- mutation_list_reports |>
  distinct() |>
  arrange(AA_pos_Ecoli, AA_mutation)

# 3.screen all rpsL sequences for existing and possible mutations:
raw_output <- screen_target_sequences(rpsL_target_sequences, rpsL_reference_Ecoli, 
                                      mutation_list, target_gene="rpsL", n_workers=6)

#save error messages:
saveRDS(raw_output[!sapply(raw_output, is.data.frame)], "./output/raw_output_errors.rds")

#save results:
raw_output <- do.call(rbind, raw_output[sapply(raw_output, is.data.frame)])
write_csv(raw_output, file = "./output/raw_output.csv")

#empty working environment to keep everything clean:
rm.all.but("globsets")

########################################################################
### Step 4: Processing and filtering of raw output                   ###
########################################################################

# input for this step:  table of reported mutations ("./output/muts.csv")
#                       screened results for all gene sequences ("./output/raw_output.csv")
#                       information of downloaded genomes from NCBI ("./output/summaries.rds)
# output for this step: filtered results for reliable sequences("./output/filtered_output.csv")
#                       summary of target sequences ("./results/summary_target_sequences.txt")
#                       plots of target sequence statistics ("target_sequence_stats_hist.pdf" & "target_sequence_stats_pairs.pdf")

# load required data:
muts <- read_csv("./output/muts.csv", show_col_types = FALSE)
raw_output <- read_csv("./output/raw_output.csv", show_col_types = FALSE)
genome_summaries <- read_rds("./output/summaries.rds")

# filter for mutations to be included in analyses:
mutation_list_reports <- filter_mutations(muts,
                                          min_n_species = globsets$min_n_species, 
                                          min_n_studies = globsets$min_n_studies)

# 2. processing and filtering raw output:
filtered_output <- raw_output |>
  process_output() |>
  # retain only data for reliable gene sequences:
  filter_output(min_seq_length = globsets$min_seq_length,
                min_alig_score = globsets$min_alig_score,
                max_core_dist = globsets$max_core_dist) |>
  # retain only data for mutations of interest:
  semi_join(mutation_list_reports, by = join_by(AA_pos_Ecoli, AA_mutation))

write_csv(filtered_output,"./output/filtered_output.csv")

# 3. analysis of extracted gene sequences and filtering:
summarise_target_sequences(genome_summaries, 
                           raw_output, 
                           filtered_output, 
                           min_seq_length = globsets$min_seq_length,
                           min_alig_score = globsets$min_alig_score,
                           max_core_dist = globsets$max_core_dist,
                           target_gene = 'rpsL',
                           file_name = "./results/summary_target_sequences.txt")
plot_target_sequences_stats(raw_output, 
                            filtered_output,
                            min_seq_length = globsets$min_seq_length,
                            min_alig_score = globsets$min_alig_score,
                            max_core_dist = globsets$max_core_dist,
                            file_names = c("./plots/target_sequence_stats_hist.pdf", "./plots/target_sequence_stats_pairs.pdf"))

#empty working environment to keep everything clean:
rm.all.but("globsets")

########################################################################
### Step 5: Analysis of results  ###
########################################################################

# input for this step:  extracted target gene sequences (./"output/rpsL_target_sequences.fa")
#                       E. coli gene reference sequence ("./data/rpsL_references.fasta")
#                       filtered output of mutation screen ("./output/filtered_output.csv")
#                       bacterial taxonomic information ("./output/bacterial_taxonomy.csv")
#                       
# output for this step: filtered results for reliable sequences("./output/filtered_output.csv")
#                       summary of screening result ("./output/summary_screen_mutations.txt")
#                       plot of predictions for each mutation across species ("./plots/mutation_screening.pdf")
#                       plot of predictions for different classes and genera ("./plots/classes_genera.pdf")
#                       table of statistics for species with multiple gene copies ("./output/multiseq_stats.csv")
#                       plot showing the statistics for species with multiple gene copies ("./plots/multicopy_stats.pdf")

# 1. load required data:
rpsL_target_sequences <- readDNAStringSet("./output/rpsL_target_sequences.fa")
rpsL_reference_Ecoli <- readDNAStringSet("./data/rpsL_references.fasta")[["rpsL_Escherichia_coli_MG1655"]]
filtered_output <- read_csv("./output/filtered_output.csv", show_col_types = FALSE)
bacterial_taxonomy <- read_csv("./data/NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI

#2. analysis of mutant screen:
plot_mutation_screen(filtered_output, file_name = "./plots/mutation_screen.pdf")
plot_classes_genera(filtered_output, bacterial_taxonomy, file_name= "./plots/classes_genera.pdf")
plot_evolvability_by_class(filtered_output, bacterial_taxonomy, file_name = "./plots/evolvability_by_class.pdf")
summarise_mutation_screen(filtered_output, target_gene = "rpsL", file_name = "./results/summary_mutation_screen.txt")
get_resistance_taxonomy(filtered_output, bacterial_taxonomy, file_path = "./output/")
make_table_intrinsic_resistance(filtered_output, file_name = "./results/predicted_resistance.csv")

#3. analyse species with multiple gene copies:
multiseq_stats <- compare_gene_copies(filtered_output, rpsL_target_sequences, rpsL_reference_Ecoli)
write_csv(multiseq_stats, "./output/multiseq_stats.csv")
#multiseq_stats <- read_csv("./output/multiseq_stats.csv", show_col_types = FALSE)
plot_multiseq_stats(multiseq_stats, "./plots/multiseq.pdf")

#empty working environment to keep everything clean:
rm.all.but("globsets")

#########################################################################
### Step 6: Phylogenetic distribution of resistance and evolvability  ###
#########################################################################

# input for this step:  filtered results for reliable sequences("./output/filtered_output.csv")
#                       original bacterial phylogenetic tree of life ("./data/bac120.nwk") and 
#                       its metadata ("./data/bac120_metadata_r214.tsv")
#                       bacterial taxonomic information ("./output/bacterial_taxonomy.csv")
#                       names of outlier species in the tree ("./output/outliers.csv")

# output for this step: subtree of original tree with tip_labels table ("./output/subtree.RData")
#                       subtree of original tree (nwk file: "./output/subtree.nwk")
#                       plot of subtree

# 1.load required files:
filtered_output <- read_csv("./output/filtered_output.csv", show_col_types = FALSE)
original_tree <- read.tree("./data/bac120.nwk") #GTDB bacterial tree of life
bacterial_taxonomy <- read_csv("./data/NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI
meta_data <- read_tsv("./data/bac120_metadata_r214.tsv", show_col_types = FALSE) #GTDB information on included species
outliers <- read_csv("./data/outliers.csv", show_col_types = FALSE) #misplaced species need to be removed later

# 2. get species-level summary of mutation screen data:
species_output <- get_species_output(filtered_output)

# 2.subset the tree based on species accessions and names:
subtree <- get_subtree(filtered_output, original_tree, meta_data, outliers)
write.tree(subtree$tree, file = "./output/subtree.nwk") 

# 3. subtree visualization:
subtree <- read.tree("./output/subtree.nwk")
# big tree of all species:
plot_subtree(subtree, species_output, bacterial_taxonomy, file_name = "./plots/phylogenies/whole_genome_tree.svg")
# smaller trees of individual clades:
plot_subtree_clades(subtree, species_output, bacterial_taxonomy, 
                    genera = c("Nocardia", "Streptomyces"),
                    families = c("Oscillospiraceae", "Francisellaceae"),
                    orders = c("Enterobacterales", "Thiotrichales"),
                    file_path = "./plots/phylogenies/")

summarise_phylogenetics(subtree, species_output, sample_n = globsets$phylo_stats_sample_n, "./results/summary_phylogenetics.txt")

#empty working environment to keep everything clean:
rm.all.but("globsets")


########################################################################
### Step 7: Codon networks                                           ###
########################################################################

# input for this step:  extracted rpoB target sequences (./"output/rpoB_target_sequences.fa")
#                       filtered output from mutation screen ("./output/filtered_output.csv")
#                       information of downloaded genomes from NCBI ("./output/summaries.rds)
#                       subtree of original tree (nwk file: "./output/subtree.nwk")
#                       
# output for this step: filtered results for reliable sequences("./output/filtered_output.csv")
#                       summary of screening result ("./output/summary_screen_mutations.txt")
#                       plots of target sequence statistics ("target_sequence_stats_hist.pdf" & "target_sequence_stats_pairs.pdf")
#                       plot of predictions for each mutation across species ("./plots/mutation_screening.pdf")
#                       plot of predictions for different genera, classes and species ("./plots/mutations_by_species.pdf")

muts <- read.csv("./output/muts.csv")
filtered_output <- read.csv("./output/filtered_output.csv")

mutation_list_reports <- filter_mutations(muts,
                                          min_n_species = globsets$min_n_species, 
                                          min_n_studies = globsets$min_n_studies,
                                          origin = "ONEMUT")

networks <- tribble(
  ~type,    ~site_order, ~pos, ~focal_codon,
  "type64", c(3, 2, 1),  516,  "GAC",
  "type64", c(3, 2, 1),  529,  "CGT",
  "type28", c(3, 2, 1),  531,  "TCG",
  "type64", c(3, 2, 1),  531,  "TCG",
  "type10", c(1, 2, 3),  522,  "TCT",
  "type28", c(2, 3, 1),  522,  "TCT",
  "type64", c(2, 3, 1),  522,  "TCT",
  "type10", c(2, 3, 1),  526,  "CAC",
  "type64", c(1, 3, 2),  513,  "CAG")

for(i in 1:nrow(networks)) {
  plot_codon_network(type = networks$type[i],
                     site_order = networks$site_order[[i]],
                     pos = networks$pos[i],
                     focal_codon = networks$focal_codon[i],
                     mutations_list = mutation_list_reports,
                     output = filtered_output,
                     lines_occupied_codons = TRUE,
                     file_path = "./plots/codon_networks/")
}

#empty working environment to keep everything clean:
rm.all.but("globsets")

########################################################################
### Step 8: Analyse mutations from Yang et al.                       ###
########################################################################

# input for this step:  table of reported mutations ("./output/muts.csv")
#                       extracted rpoB target sequences (./"output/rpoB_target_sequences.fa")
#                       screened results for all rpoB sequences ("./output/raw_output.csv")
#                       reported frequencies of rpoB mutations under RIF, 
#                           from Yang et al. (2023) ("./data/Yang2023_Fig3a_data_RIF_log2FC.csv")
#                       subtree of original tree (nwk file: "./output/subtree.nwk")
#                       
# output for this step: filtered results for reliable sequences("./output/filtered_output_YangOnly.csv")
#                       summary of screening result ("./output/summary_screen_mutations.txt")
#                       plots of target sequence statistics ("target_sequence_stats_hist.pdf" & "target_sequence_stats_pairs.pdf")
#                       plot of predictions for each mutation across species ("./plots/mutation_screening.pdf")
#                       plot of predictions for different classes and genera ("./plots/mutations_by_species.pdf")

# load required data:
muts <- read_csv("./output/muts.csv", show_col_types = FALSE)
rpoB_target_sequences <- readDNAStringSet("./output/rpoB_target_sequences.fa")
raw_output <- read_csv("./output/raw_output.csv", show_col_types = FALSE)
subtree <- read.tree("./output/subtree.nwk")
bacterial_taxonomy <- read_csv("./data/NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI
mutation_list_Yang <- import_Yang_muts(file_name_RIF_log2FC = "./data/Yang2023_Fig3a_data_RIF_log2FC.csv",
                                       min_log2_fc = globsets$Yang_min_log2_fc)

# 1. filtering final output for reliable sequences
filtered_output <- raw_output |>
  process_output() |>
  # retain only data for reliable rpoB sequences:
  filter_output(min_seq_length = globsets$min_seq_length,
                min_alig_score = globsets$min_alig_score,
                max_core_dist = globsets$max_core_dist)

# 2.filter to only retain mutations that have not been reported but are in Yang et al.:
mutation_list_reports <- filter_mutations(muts,
                                          min_n_species = globsets$min_n_species, 
                                          min_n_studies = globsets$min_n_studies,
                                          origin = "ONEMUT")

filtered_output <- filtered_output |>
  anti_join(mutation_list_reports, by = join_by(AA_pos_Ecoli, AA_mutation))

write_csv(filtered_output,"./output/filtered_output_YangOnly.csv")
#filtered_output <- read_csv("./output/filtered_output_YangOnly.csv", show_col_types = FALSE)

#3. analysis of filtered results:
plot_reported_vs_Yang(mutation_list_reports, mutation_list_Yang, file_name = "./plots/YangOnly/reported_vs_Yang_mutations.pdf")
plot_mutation_screen(filtered_output, file_name = "./plots/YangOnly/mutation_screening_YangOnly.pdf")
plot_classes_genera(filtered_output, bacterial_taxonomy, file_name= "./plots/YangOnly/mutations_by_species_YangOnly.pdf")
summarise_mutation_screen(filtered_output,
                          file_name = "./results/summary_mutation_screen_YangOnly.txt",
                          subtitle = "(additional mutations from Yang et al. (2023))")

#4. phylogenetic tree:
species_output <- get_species_output(filtered_output)
plot_subtree(subtree, species_output, bacterial_taxonomy, file_name = "./plots/YangOnly/whole_genome_tree_YangOnly.pdf")

#empty working environment to keep everything clean:
rm.all.but("globsets")


########################################################################
### Step 9: Analyse all mutations (reports & from Yang et al.        ###
########################################################################

# input for this step:  extracted rpoB target sequences (./"output/rpoB_target_sequences.fa")
#                       screened results for all rpoB sequences ("./output/raw_output.csv")
#                       subtree of original tree (nwk file: "./output/subtree.nwk")
#                       
# output for this step: summary of screening result ("./output/summary_screen_mutations.txt")
#                       plot of predictions for each mutation across species ("./plots/mutation_screening_AllMutations.pdf")
#                       plot of predictions for different classes and genera ("./plots/mutations_by_species_AllMutations.pdf")

# load required data:
rpoB_target_sequences <- readDNAStringSet("./output/rpoB_target_sequences.fa")
raw_output <- read_csv("./output/raw_output.csv", show_col_types = FALSE)
subtree <- read.tree("./output/subtree.nwk")
bacterial_taxonomy <- read_csv("./data/NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI

# 1. filtering final output for reliable sequences
filtered_output <- raw_output |>
  process_output() |>
  # retain only data for reliable rpoB sequences:
  filter_output(min_seq_length = globsets$min_seq_length,
                min_alig_score = globsets$min_alig_score,
                max_core_dist = globsets$max_core_dist)

#3. analysis of filtered results:
plot_mutation_screen(filtered_output, file_name = "./plots/AllMutations/mutation_screening_AllMutations.pdf")
plot_classes_genera(filtered_output, bacterial_taxonomy, file_name= "./plots/AllMutations/mutations_by_species_AllMutations.pdf")
summarise_mutation_screen(filtered_output,
                          file_name = "./results/summary_mutation_screen_AllMutations.txt",
                          subtitle = "(reported mutations and mutations from Yang et al. (2023) pooled)")

#4. phylogenetic tree:
species_output <- get_species_output(filtered_output)
plot_subtree(subtree, species_output, bacterial_taxonomy, file_name = "./plots/AllMutations/whole_genome_tree_AllMutations.pdf")

#empty working environment to keep everything clean:
rm.all.but("globsets")

########################################################################
### Step 10: Plotting rpoB structure and sequence conservation       ###
########################################################################

# load required data:
rpoB_reference_Ecoli <- readDNAStringSet("./data/rpoB_references.fasta")[["rpoB_Escherichia_coli_MG1655"]]
filtered_output <- read_csv("./output/filtered_output.csv", show_col_types = FALSE)
filtered_targets <- filtered_output |> pull(target_name) |> unique()
rpoB_target_sequences <- readDNAStringSet("./output/rpoB_target_sequences.fa")[filtered_targets]

mutations <- read_csv("./output/muts.csv", show_col_types = FALSE) |>
  filter_mutations(min_n_species = globsets$min_n_species, 
                   min_n_studies = globsets$min_n_studies)

# calculate conservation/diversity scores along the rpoB sequence:
set.seed(globsets$random_seed)
cons <- get_conservation(rpoB_target_sequences, rpoB_reference_Ecoli, n_rnd = 1e5)
save(cons, file = "./output/cons.RData")
summarise_conservation(cons,
                       file_name = "./results/summary_conservation.txt")
plot_cons(cons, 
          pos = mutations |> pull(AA_pos_Ecoli) |> unique(),
          dist_type = "hamming",
          file_name = "./plots/AA_conservation_hamming.pdf")
plot_cons(cons, 
          pos = mutations |> pull(AA_pos_Ecoli) |> unique(),
          dist_type = "grantham",
          file_name = "./plots/AA_conservation_grantham.pdf")

# download PDB file for E. coli rpoB:
pdb <- read.pdb("5uac")
# Select beta chain (all chains are duplicated in crystal structure):
chain_selection <- atom.select(pdb, chain = "C")
# Extract the selected chain:
pdb_chain <- trim.pdb(pdb, chain_selection)
# Save the selected chain to a PDB file
write.pdb(pdb_chain, file = "./data/5uac_chainC.pdb")
# Produce html file of protein structure with distance to E. coli indicated:
visualise_RpoB_structure(file_pdb = "./data/5uac_chainC.pdb",
                         pos = mutations |> pull(AA_pos_Ecoli) |> unique(),
                         mut_colour_variable = cons$means$grantham_Ecoli,
                         file_html = "./plots/rpoB_structure.html")
quarto_render("plots/rpoB_structure_embedding.qmd")

########################################################################
### Step 11: Compile final report                                    ###
########################################################################

# collating all individual summary files and rendering them as a single pdf file:
render_summary(summaries = c("reported_mutations", 
                             "target_sequences", 
                             "mutation_screen", 
                             "phylogenetics",
                             "mutation_screen_YangOnly",
                             "mutation_screen_AllMutations",
                             "conservation"),
              preamble = "./data/summary_preamble.qmd",
              summaries_path = "./results")

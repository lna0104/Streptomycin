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
library(rentrez)

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
library(wesanderson)
library(ggrepel)
library(patchwork)
library(cowplot) 
library(ggraph)
library(tidygraph)
library(ggh4x)
library(igraph)



source("R/util.R")
source("R/bioinformatics.R")
source("R/analyses.R")
source("R/plotting.R")
source("R/reports.R")
source("R/phylogenetics.R")
source("R/codon_networks.R")
source("R/structure.R")

# Global settings and parameters:
globsets <- list(
  min_n_studies = 3, # minimum number of studies that a mutation needs to be reported in for inclusion
  min_n_species = 3, # minimum number of species that a mutation needs to be reported in for inclusion
  min_seq_length = 300, # minimum length of included gene target sequences
  min_alig_score = -Inf, # minimum alignment score (with E. coli) of included gene target sequences
  max_core_dist = 40, # maximum Levenshtein distance between E. coli core gene region to corresponding target region
  phylo_stats_sample_n = 5000, # number of species to sample for phylogenetics statistics
  random_seed = 22)
options(nwarnings = 10000)

set.seed(globsets$random_seed)

########################################################################
### Step 1: Processing reference files                              ###
########################################################################
# Read the mutations data from a CSV file
muts<-read.csv("./data/reported_mutations.csv") |>
  filter(Gene=="rpsL")

# Select unique species from the mutations data and create a gene reference data frame
selected_muts <- muts |>
  select(Species) |>  # Select the 'Species' column
  unique() |>  # Remove duplicate species names
  mutate(
    ID = NA,  # Initialize the 'ID' column with NA
    Strain = NA,  # Initialize the 'Strain' column with NA
    RefStrain = NA,  # Initialize the 'RefStrain' column with NA
    Gene = "rpsL",  # Add the gene name as 'rpsL'
    Assembly_ID = NA,  # Initialize the 'Assembly_ID' column with NA
    NCBI_Reference_sequence = NA  # Initialize the 'NCBI_Reference_sequence' column with NA
  )

summaries<-list()

# Loop through each species in the selected_muts data frame
for (i in 1:nrow(selected_muts)){
  species = selected_muts$Species[i]
  cat(paste0("Working on i = ",i, " - ", species, " species\n"))
  
  {
    # Search for reference genomes for the species in the NCBI assembly database
    # Check if the species is Escherichia coli and include strain in the search term
    if (grepl("Escherichia coli", species, ignore.case = TRUE)) {
      term_ser <- sprintf('"%s"[ORGN] AND "MG1655"[Strain] AND ("latest refseq"[filter] AND "reference genome"[filter])', species)
    } else {
      term_ser <- sprintf('"%s"[ORGN] AND ("latest refseq"[filter] AND "reference genome"[filter])', species)
    }
    assembly_search_results <- get_summaries(db="assembly", term=term_ser)
    # Check if any assembly search results were found
    if (length(assembly_search_results)!=0){
      # If found, select the first assembly result as the reference genome
      assembly_search_final_result<-assembly_search_results[[1]]
      selected_muts$RefStrain[i]<-TRUE
    }else{
      # If no reference genome was found, search for the species in general
      term_ser<- sprintf('"%s"[ORGN] AND "latest refseq"[filter]', species)
      assembly_search_results <- get_summaries(db="assembly", term=term_ser)
      if (length(assembly_search_results)==0){
        cat(paste0("Cannot find assembly data for ", species, ". Try search species name in NCBI."))
        next
      }
      assembly_search_final_result<-assembly_search_results[[1]]
      selected_muts$RefStrain[i]<-FALSE
    }
    # Store the assembly summary in the 'summaries' list using the assembly UID as the key
    summaries[[assembly_search_final_result$uid]]<-assembly_search_final_result
    selected_muts$ID[i]<-assembly_search_final_result$uid
  }
  # Store the assembly name in the 'Assembly_ID' column
  selected_muts$Assembly_ID[i]<-assembly_search_final_result$assemblyname
  
  {
    # Fetch links to nucleotide sequence data for the selected assembly
    link_to_nuc <- entrez_link(dbfrom = "assembly", db = "nuccore", id = assembly_search_final_result$uid)
    n_link_ref <- length(link_to_nuc$links$assembly_nuccore_refseq)
    n_link <- length(link_to_nuc$links$assembly_nuccore)
    
    # Check if any nucleotide sequences are found
    if (n_link == 0) {
      cat("No nucleotide information found! Skipping\n")
      next
    }
    
    # Function to extract strain and reference sequence
    extract_strain_info <- function(nuc_result) {
      strain <- nuc_result$strain
      ref_seq <- nuc_result$accessionversion
      
      # If the species is E. coli, force strain to be MG1655
      if (grepl("Escherichia coli", species, ignore.case = TRUE)) {
        strain <- "MG1655"
      }
      
      return(list(strain = strain, ref_seq = ref_seq))
    }
    
    # If no refseq sequences were found
    if (n_link_ref == 0) {
      nuc_search_results <- entrez_summary(db = "nuccore", id = link_to_nuc$links$assembly_nuccore)
      
      if (n_link > 1) {
        strain_info <- extract_strain_info(nuc_search_results[[1]])
      } else if (n_link == 1) {
        strain_info <- extract_strain_info(nuc_search_results)
      }
    }
    
    # If refseq sequences are found
    if (n_link_ref >= 1) {
      nuc_search_results <- entrez_summary(db = "nuccore", id = link_to_nuc$links$assembly_nuccore_refseq)
      
      if (n_link_ref > 1) {
        strain_info <- extract_strain_info(nuc_search_results[[1]])
      } else {
        strain_info <- extract_strain_info(nuc_search_results)
      }
    }
    
    # Assign extracted strain and reference sequence
    selected_muts$Strain[i] <- strain_info$strain
    selected_muts$NCBI_Reference_sequence[i] <- strain_info$ref_seq
  }
}

#Generate the fasta name
added_name_muts <- selected_muts |>
  mutate(Species_name = gsub(" ", "_", Species),  # Replace spaces with underscores
         Strain_name = gsub(" ", "_", Strain),
         FASTA_name = paste0(Gene, "_", Species_name, "_", Strain_name)) |>
  select(-Species_name, Strain_name) 

# Save the table as a CSV file
write.csv(added_name_muts, "data/rpsL_references.csv", row.names = FALSE)

# Initialize an empty list to store the sequences from all downloads
all_sequences <- list()
# Loop through each genome summary in the 'summaries' list
for (j in 1:length(summaries)) {
  id <- names(summaries[j])
  cat(paste0("Downloading genome for species #", j, "/", length(summaries), "\n..."))
  file_path <- download_file(summaries[[j]], dir = "output/rpsL_references")
  fasta_file <- sub("\\.gz$", "", file_path)
  sequences <- readDNAStringSet(fasta_file)
  matching_seq <- sequences[grep("rpsL|30S ribosomal protein S12($|\\])", names(sequences))]
  if (length(matching_seq) == 0){
    cat(paste0("No rpsL found for species ID ", id, "\n"))
    next
  } else if (length(matching_seq) > 1) {
    cat("There are ", length(matching_seq), "copies of rpsL. \n")
    matching_seq <- matching_seq[1]
  }
  names(matching_seq) <- added_name_muts[added_name_muts$ID == id, ]$FASTA_name
  all_sequences <- c(all_sequences, matching_seq)
}
# Combine all sequences into a single DNAStringSet object
combined_sequences <- do.call(c, all_sequences)

# Write the combined sequences to a new FASTA file
writeXStringSet(combined_sequences, "data/rpsL_references.fasta")

#empty working environment to keep everything clean
rm.all.but("globsets")

########################################################################
### Step 2: Processing table of known STR resistance mutations       ###
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

compareMutationsToRef<-function(seqs, muts){
  muts<-muts |>
    mutate(Warning2 = NA)
  for (k in 1:nrow(muts)){
    seq_name<-muts$RefSeq_ID[k]
    AA_seq<-translate(seqs[seq_name])[[1]]
    AA_pos<-muts$AA_pos[k]
    if (as.character(AA_seq[AA_pos]) != muts$AA_original[k] & as.character(AA_seq[AA_pos]) != muts$AA_mutation[k]){
      muts$Warning2[k]<-"AA_original or AA_mutation inconsistent with AA_pos"
    }
  }
  return(muts)
}

added_warnings_muts<-compareMutationsToRef(seqs, muts)

write_csv(added_warnings_muts, "./output/muts.csv")
#3 summary and plot of reported mutations
plot_reported_mutations(added_warnings_muts, file_name = "./plots/reported_mutations_original.pdf", n_frequency = 3) # returns frequent reported mutations, positions and species
summarise_reported_mutations(added_warnings_muts, file_name = "./results/summary_reported_mutations_original.txt") # returns a text message summarizing previous reports

#Manually check all warnings and correct mutations
checked_muts<-read.csv("./output/checked_muts.csv")

#3 summary and plot of reported mutations
plot_reported_mutations(checked_muts, file_name = "./plots/reported_mutations_manualfix.pdf", n_frequency = 3) # returns frequent reported mutations, positions and species
summarise_reported_mutations(checked_muts, file_name = "./results/summary_reported_mutations_manualfix.txt") # returns a text message summarizing previous reports

#empty working environment to keep everything clean
rm.all.but("globsets")

##############################################################################
### Step 3: Retrieving target sequences from all bacterial reference genomes ###
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
                                              target_protein = "30S ribosomal subunit protein S12|30S ribosomal protein S12($|[^'])")
writeXStringSet(DNAStringSet(rpsL_target_sequences), filepath = "output/rpsL_target_sequences.fa")

# download and extract taxonomy information for downloaded genomes
download_taxonomy(summaries, output_file = "./data/NCBI_taxonomy.csv")

#empty working environment to keep everything clean
rm.all.but("globsets")

#########################################################################
### Step 4: Checking all target sequences for reported mutations      ###
#########################################################################

# input for this step:  processed table of reported mutations after manually modify all warnings ("./output/checked_muts.csv")
#                       reported frequencies of rpsL mutations under STR 
#                       E. coli rspL reference sequence (from "./data/rspL_references.fasta")
#                       extracted rspL target sequences ("./output/rspL_target_sequences.fa")
# output for this step: Table of screening results, including existing and evolvable mutations 
#                           for all target sequences ("./output/raw_output.csv")

# 1.load required data:
muts <- read.csv("./output/checked_muts.csv") 
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
### Step 5: Processing and filtering of raw output                   ###
########################################################################

# input for this step:  table of reported mutations ("./output/checked_muts.csv")
#                       screened results for all gene sequences ("./output/raw_output.csv")
#                       information of downloaded genomes from NCBI ("./output/summaries.rds)
# output for this step: filtered results for reliable sequences("./output/filtered_output.csv")
#                       summary of target sequences ("./results/summary_target_sequences.txt")
#                       plots of target sequence statistics ("target_sequence_stats_hist.pdf" & "target_sequence_stats_pairs.pdf")

# load required data:
muts <- read_csv("./output/checked_muts.csv", show_col_types = FALSE) 
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
### Step 6: Analysis of results  ###
########################################################################

# input for this step:  extracted target gene sequences (./"output/rpsL_target_sequences.fa")
#                       E. coli gene reference sequence ("./data/rpsL_references.fasta")
#                       filtered output of mutation screen ("./output/filtered_output.csv")
#                       GTDB bacterial taxonomic information ("./output/gtdb_taxonomy.csv")
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
# bacterial_taxonomy <- read_csv("./data/NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI
# meta_data <- read_tsv("./data/bac120_metadata.tsv", show_col_types = FALSE) #GTDB information on included species
# meta_data_parsed <- meta_data |>
#   select(gtdb_taxonomy) |>
#   separate(gtdb_taxonomy, 
#            into = c("domain", "phylum", "class", "order", "family", "genus", "species"), 
#            sep = ";", 
#            remove = FALSE) %>%
#   mutate(across(domain:species, ~ sub("^[a-z]__","", .))) |>
#   select(genus, family, order, class, phylum) |>
#   distinct()
# write_csv(meta_data_parsed, "./data/gtdb_taxonomy.csv") 
# Generate gtdb_taxonomy file from GTDB metadata
gtdb_taxonomy <- read_csv("./data/gtdb_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from gtdb
genus_variants <- read_csv("./output/variants_gtdb_taxonomy.csv", show_col_types = FALSE)


#2. analysis of mutant screen:
plot_mutation_screen(filtered_output, file_name = "./plots/mutation_screen.pdf")
plot_classes_genera(filtered_output, gtdb_taxonomy, genus_variants, file_name= "./plots/classes_genera.pdf")
# plot_classes_genera(filtered_output, gtdb_taxonomy, genus_variants, file_name= "./plots/classes_genera.svg")
plot_evolvability_by_class(filtered_output, gtdb_taxonomy, genus_variants, file_name = "./plots/evolvability_by_class.pdf")
summarise_mutation_screen(filtered_output, target_gene = "rpsL", file_name = "./results/summary_mutation_screen.txt")
get_resistance_taxonomy(filtered_output, gtdb_taxonomy, genus_variants,  file_path = "./output/")
make_table_intrinsic_resistance(filtered_output, file_name = "./results/predicted_resistance.csv")

#3. analyse species with multiple gene copies:
multiseq_stats <- compare_gene_copies(filtered_output, rpsL_target_sequences, rpsL_reference_Ecoli)
write_csv(multiseq_stats, "./output/multiseq_stats.csv")
#multiseq_stats <- read_csv("./output/multiseq_stats.csv", show_col_types = FALSE)
plot_multiseq_stats(multiseq_stats, "./plots/multiseq.pdf")

#empty working environment to keep everything clean:
rm.all.but("globsets")

#########################################################################
### Step 7: Phylogenetic distribution of resistance and evolvability  ###
#########################################################################

# input for this step:  filtered results for reliable sequences("./output/filtered_output.csv")
#                       original bacterial phylogenetic tree of life ("./data/bac120.nwk") and 
#                       its metadata ("./data/bac120_metadata.tsv")
#                       GTDB bacterial taxonomic information ("./output/gtdb_taxonomy.csv")

# output for this step: subtree of original tree with tip_labels table ("./output/subtree.RData")
#                       subtree of original tree (nwk file: "./output/subtree.nwk")
#                       plot of subtree

# 1.load required files:
filtered_output <- read_csv("./output/filtered_output.csv", show_col_types = FALSE)
original_tree <- read.tree("./data/bac120.nwk") #GTDB bacterial tree of life
# original_tree <- read.tree("./data/bac120.tree") #GTDB bacterial tree of life
# bacterial_taxonomy <- read_csv("./data/NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI
meta_data <- read_tsv("./data/bac120_metadata.tsv", show_col_types = FALSE) #GTDB information on included species
gtdb_taxonomy <- read_csv("./data/gtdb_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from gtdb
# outliers <- if (file.exists("./data/outliers.csv")) {
#   read_csv("./data/outliers.csv", show_col_types = FALSE)
# } else {
#   NULL
# }
genus_variants <- read_csv("./output/variants_gtdb_taxonomy.csv", show_col_types = FALSE)

# 2. get species-level summary of mutation screen data:
species_output <- get_species_output(filtered_output)

# 2.subset the tree based on species accessions and names:
subtree <- get_subtree(filtered_output, original_tree, meta_data)
write.tree(subtree$tree, file = "./output/subtree.nwk") 

# 3. subtree visualization:
subtree <- read.tree("./output/subtree.nwk")
# big tree of all species:
plot_subtree(subtree, species_output, gtdb_taxonomy, genus_variants, file_name = "./plots/phylogenies/whole_genome_tree.svg")
# smaller trees of individual clades:
plot_subtree_clades(subtree, species_output, gtdb_taxonomy, genus_variants, 
                    genera = c("Sphingomonas"),
                    families = c("Devosiaceae", "Mycobacteriaceae"),
                    orders = c("Pirellulales", "Sphingomonadales", "Rickettsiales"),
                    classes = c("Planctomycetia", "Alphaproteobacteria","Coriobacteriia"),
                    file_path = "./plots/phylogenies/")

summarise_phylogenetics(subtree, species_output, sample_n = globsets$phylo_stats_sample_n, "./results/summary_phylogenetics.txt")

#empty working environment to keep everything clean:
rm.all.but("globsets")


########################################################################
### Step 8: Codon networks                                         ###
########################################################################

# input for this step:  extracted target sequences (./"output/rpsL_target_sequences.fa")
#                       filtered output from mutation screen ("./output/filtered_output.csv")
#                       information of downloaded genomes from NCBI ("./output/summaries.rds)
#                       subtree of original tree (nwk file: "./output/subtree.nwk")
#                       
# output for this step: filtered results for reliable sequences("./output/filtered_output.csv")
#                       summary of screening result ("./output/summary_screen_mutations.txt")
#                       plots of target sequence statistics ("target_sequence_stats_hist.pdf" & "target_sequence_stats_pairs.pdf")
#                       plot of predictions for each mutation across species ("./plots/mutation_screening.pdf")
#                       plot of predictions for different genera, classes and species ("./plots/mutations_by_species.pdf")

muts <- read.csv("./output/checked_muts.csv")
filtered_output <- read.csv("./output/filtered_output.csv")

mutation_list_reports <- filter_mutations(muts,
                                          min_n_species = globsets$min_n_species, 
                                          min_n_studies = globsets$min_n_studies,
                                          origin = "ONEMUT")

networks <- tribble(
  ~type,    ~site_order, ~pos, ~focal_codon,
  "type64", c(3, 2, 1),  43,  "AAG",
  "type28", c(3, 2, 1),  43,  "AAG",
  "type64", c(3, 2, 1),  86,  "CGT",
  "type28", c(3, 2, 1),  86,  "CGT",
  "type64", c(3, 2, 1),  92,  "GGT",
  "type28", c(3, 2, 1),  92,  "GGT",
  "type64", c(3, 2, 1),  88,  "AAG",
  "type28", c(3, 2, 1),  88,  "AAG",
  )

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
### Step 9: Plotting rpsL structure and sequence conservation       ###
########################################################################

# load required data:
rpsL_reference_Ecoli <- readDNAStringSet("./data/rpsL_references.fasta")[["rpsL_Escherichia_coli_MG1655"]]
filtered_output <- read_csv("./output/filtered_output.csv", show_col_types = FALSE)
filtered_targets <- filtered_output |> pull(target_name) |> unique()
rpsL_target_sequences <- readDNAStringSet("./output/rpsL_target_sequences.fa")[filtered_targets]

mutations <- read_csv("./output/checked_muts.csv", show_col_types = FALSE) |>
  filter_mutations(min_n_species = globsets$min_n_species, 
                   min_n_studies = globsets$min_n_studies)

# calculate conservation/diversity scores along the rpsL sequence:
set.seed(globsets$random_seed)
cons <- get_conservation(rpsL_target_sequences, rpsL_reference_Ecoli, n_rnd = 1e5, n_workers=10)
save(cons, file = "./output/cons.RData")
summarise_conservation(cons, 
                       target_gene = "rpsL",
                       file_name = "./results/summary_conservation.txt")
plot_cons(cons, 
          pos = mutations |> pull(AA_pos_Ecoli) |> unique(),
          pos_range = c(40,100),
          n_plots = 1,
          dist_type = "hamming",
          file_name = "./plots/AA_conservation_hamming.pdf")
plot_cons(cons, 
          pos = mutations |> pull(AA_pos_Ecoli) |> unique(),
          pos_range = c(40,100),
          n_plots = 1,
          dist_type = "grantham",
          file_name = "./plots/AA_conservation_grantham.pdf")

# download PDB file for E. coli rpsL:
pdb <- read.pdb("8cgj")
# pdb <- read.pdb("8cai")
# Select L chain (Small ribosomal subunit protein uS12):
# Select A chain (16S rRNA where streptomycin bind to)
chain_selection <- atom.select(pdb, chain = c("L", "A"))
# Extract the selected structure:
pdb_chain <- trim.pdb(pdb, chain_selection)
# Save the selected structure to a PDB file
write.pdb(pdb_chain, file = "./data/8cgj_chain.pdb")
# Produce html file of protein structure with distance to E. coli indicated:
# load("./output/cons.RData")
visualise_rpsL_structure(file_pdb = "./data/8cgj_chain.pdb",
                         chain_ids = c("L","A"),
                         ligand_resid = "5I0",
                         pos = mutations |> pull(AA_pos_Ecoli) |> unique(),
                         mut_colour_variable = cons$means$grantham_Ecoli,
                         file_html = "./plots/rpsL_structure.html")
quarto_render("plots/rpsL_structure_embedding.qmd")

########################################################################
### Step 10: Compile final report                                    ###
########################################################################

# collating all individual summary files and rendering them as a single pdf file:
render_summary(summaries = c("reported_mutations_manualfix", 
                             "target_sequences", 
                             "mutation_screen", 
                             "phylogenetics",
                             "conservation"),
              preamble = "./data/summary_preamble.qmd",
              summaries_path = "./results")

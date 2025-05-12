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



source("R/util.R")
source("R/bioinformatics_nt.R")
source("R/bioinformatics.R")
source("R/analyses.R")
source("R/plotting_nt.R")
source("R/plotting.R")
source("R/reports_nt.R")
source("R/reports.R")
source("R/phylogenetics.R")
source("R/codon_networks.R")
source("R/structure.R")
source("R/process_nt_data.R")

# Global settings and parameters:
globsets <- list(
  min_n_studies = 3, # minimum number of studies that a mutation needs to be reported in for inclusion
  min_n_species = 3, # minimum number of species that a mutation needs to be reported in for inclusion
  min_seq_length = 1200, # minimum length of included gene target sequences
  min_alig_score = -Inf, # minimum alignment score (with E. coli) of included gene target sequences
  max_core_dist = 120, # maximum Levenshtein distance between E. coli core gene region to corresponding target region
  phylo_stats_sample_n = 5000, # number of species to sample for phylogenetics statistics
  random_seed = 22)
options(nwarnings = 10000)

set.seed(globsets$random_seed)

########################################################################
### Step 1: Processing reference files                              ###
########################################################################
# Read the mutations data from a CSV file
muts<-read.csv("./data/reported_mutations.csv") |>
  filter(Gene=="rrs")

# Select unique species from the mutations data and create a gene reference data frame
selected_muts <- muts |>
  select(Species) |>  # Select the 'Species' column
  unique() |>  # Remove duplicate species names
  mutate(
    ID = NA,  # Initialize the 'ID' column with NA
    Strain = NA,  # Initialize the 'Strain' column with NA
    RefStrain = NA,  # Initialize the 'RefStrain' column with NA
    Gene = "rrs",  # Add the gene name as 'rpsL'
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
write.csv(added_name_muts, "data/rrs_references.csv", row.names = FALSE)

# Initialize an empty list to store the sequences from all downloads
all_sequences <- list()
# Loop through each genome summary in the 'summaries' list
for (j in 1:length(summaries)) {
  id <- names(summaries[j])
  cat(paste0("Downloading genome for species #", j, "/", length(summaries), "\n..."))
  file_path <- download_file(summaries[[j]], suffix = "_rna_from_genomic.fna.gz", dir = "output/rrs_references")
  fasta_file <- sub("\\.gz$", "", file_path)
  sequences <- readDNAStringSet(fasta_file)
  matching_seq <- sequences[grep("rrs|16S ribosomal RNA($|\\])", names(sequences))]
  if (length(matching_seq) == 0){
    cat(paste0("No rrs found for species ID ", summaries[[j]]$speciesname, "\n"))
    next
  } else if (length(matching_seq) > 1) {
    cat("There are ", length(matching_seq), "copies of rpsL of", summaries[[j]]$speciesname, "\n")
    #distance_dist <- stringDist(matching_seq, method = "levenshtein")
    least_different_name <- least_different_sequence(matching_seq)
    matching_seq <- matching_seq[least_different_name]
  }
  names(matching_seq) <- added_name_muts[added_name_muts$ID == id, ]$FASTA_name
  all_sequences <- c(all_sequences, matching_seq)
}
# Combine all sequences into a single DNAStringSet object
combined_sequences <- do.call(c, all_sequences)

# Write the combined sequences to a new FASTA file
writeXStringSet(combined_sequences, "data/rrs_references.fasta")

#empty working environment to keep everything clean
rm.all.but("globsets")


########################################################################
### Step 2: Processing table of known STR resistance mutations       ###
########################################################################

# input for this step: table of rrs reference sequences ("./data/rrs_references.csv")
#                      fasta file with rrs reference sequences (""./data/rrs_references.fasta"")
#                      table of reported STR mutations ("./data/reported_mutations.csv")
# output for this step: coordinates ("./output/coordinates.RData"), 
#                       table of reported mutations ("./output/muts.csv")
#                       a plot showing distributions of reported mutations ("./plots/reported_mutations.pdf)
#                       a summary of the reported mutations ("./results/summary_reported_mutations.txt")

# 1 get rpsL ref sequences coordination of all reference sequences against the one in Escherichia coli

#load reference sequences and their information 
refs <- read_csv("./data/rrs_references.csv", show_col_types = FALSE)
seqs <- readDNAStringSet("./data/rrs_references.fasta")
mutation_list_reports <- read_csv("./data/reported_mutations.csv", show_col_types = FALSE)|>
  filter(Gene=="rrs")

# Check whether the above files have been changed and hence the coordinates need to be updated 
if(file.exists("./data/fastahash_rrs.Rds") && as.character(openssl::sha1(file("./data/rrs_references.fasta"))) == readRDS("./data/fastahash_rrs.Rds")){
  print("Sequences file has not changed, loading original coordinates")
  load(file = "./output/coordinates_rrs.RData")
} else {
  old_fastahash <- as.character(openssl::sha1(file("./data/rrs_references.fasta")))
  saveRDS(old_fastahash, 
       file = "./data/fastahash_rrs.Rds")
  print("Sequences file has changed, regenerating coordinates")
  #get coordinates
  coordinates <- ALJEbinf::getAllCoordinates(seqs, "rrs_Escherichia_coli_MG1655")
  save(coordinates, file = "./output/coordinates_rrs.RData")
}

#2 load and complete table of reported mutations:
muts <- mutation_list_reports |>
  fillMutationsTableNt(refs, seqs, coordinates) |>
  filter(!is.na(Nt_mut_name_Ecoli))

write_csv(muts, "./output/muts_rrs.csv")
#3 summary and plot of reported mutations
plot_reported_mutations_nt(muts, file_name = "./plots/rrs_reported_mutations.pdf", n_frequency = 2) # returns frequent reported mutations, positions and species
summarise_reported_mutations_nt(muts, file_name = "./results/rrs_summary_reported_mutations.txt") # returns a text message summarizing previous reports

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
# term_rep <- '("Escherichia coli"[Organism] OR "Mycobacterium"[Organism]) AND ("latest refseq"[filter] AND "representative genome"[filter] AND "refseq has annotation"[Properties])'#representative

term_rep <- '("Bacteria"[Organism] OR bacteria[All Fields]) AND ("latest refseq"[filter] AND "representative genome"[filter] AND "refseq has annotation"[Properties])'#representative
#search the entire database using the defined term
summaries_rep <- get_summaries(db, term_rep)

# define search term for reference genomes: Bacterial genomes at all assembly levels with annotation 
# term_ref <- '("(Escherichia coli"[Organism] OR "Mycobacterium"[Organism]) AND ("latest refseq"[filter] AND "reference genome"[filter] AND "refseq has annotation"[Properties])'

term_ref <- '("Bacteria"[Organism] OR bacteria[All Fields]) AND ("latest refseq"[filter] AND "reference genome"[filter] AND "refseq has annotation"[Properties])'
# search the entire database using the defined term
summaries_ref <- get_summaries(db, term_ref)

# combine summaries and save them
summaries <- c(summaries_ref, summaries_rep)
saveRDS(summaries, file = "./output/rrs_summaries.rds")

# download genomes
download_files(summaries, dir = "output/rrs_genomes", suffix = "_rna_from_genomic.fna.gz")

# extract target sequences and save them as rds and fasta
rrs_target_sequences <- get_target_sequences_nt(summaries, 
                                                dir = "output/rrs_genomes", 
                                                suffix = "_rna_from_genomic.fna", 
                                                target_gene = 'rrs', 
                                                target_rna = "16S ribosomal RNA")
writeXStringSet(DNAStringSet(rrs_target_sequences), filepath = "output/rrs_target_sequences.fa")

# download and extract taxonomy information for downloaded genomes
download_taxonomy(summaries, output_file = "./data/rrs_NCBI_taxonomy.csv")

#empty working environment to keep everything clean
rm.all.but("globsets")

#########################################################################
### Step 4: Checking all target sequences for reported mutations      ###
#########################################################################

# input for this step:  processed table of reported mutations after manually modify all warnings ("./output/muts_rrs.csv")
#                       reported frequencies of rrs mutations under STR 
#                       E. coli rrs reference sequence (from "./data/rrs_references.fasta")
#                       extracted rrs target sequences ("./output/rrs_target_sequences.fa")
# output for this step: Table of screening results, including existing and evolvable mutations 
#                           for all target sequences ("./output/rrs_raw_output.csv")

# 1.load required data:
muts <- read.csv("./output/muts_rrs.csv") 
rrs_target_sequences <- readDNAStringSet("./output/rrs_target_sequences.fa")
rrs_reference_Ecoli <- readDNAStringSet("./data/rrs_references.fasta")[["rrs_Escherichia_coli_MG1655"]]

# 2.make a list of reliable mutations to be screened:
mutation_list_reports <- filter_mutations_nt(muts,
                                          min_n_species = globsets$min_n_species, 
                                          min_n_studies = globsets$min_n_studies)
mutation_list <- mutation_list_reports |>
  distinct() |>
  arrange(Nt_pos_Ecoli, Nt_mutation)

# 3.screen all rrs sequences for existing and possible mutations:
raw_output <- screen_target_sequences_nt(rrs_target_sequences, rrs_reference_Ecoli, 
                                      mutation_list, target_gene="rrs", n_workers=6)

#save error messages:
saveRDS(raw_output[!sapply(raw_output, is.data.frame)], "./output/rrs_raw_output_errors.rds")

#save results:
raw_output <- do.call(rbind, raw_output[sapply(raw_output, is.data.frame)])
write_csv(raw_output, file = "./output/rrs_raw_output.csv")

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
muts <- read_csv("./output/muts_rrs.csv", show_col_types = FALSE) 
raw_output <- read_csv("./output/rrs_raw_output.csv", show_col_types = FALSE)
genome_summaries <- read_rds("./output/rrs_summaries.rds")

# filter for mutations to be included in analyses:
mutation_list_reports <- filter_mutations_nt(muts,
                                          min_n_species = globsets$min_n_species, 
                                          min_n_studies = globsets$min_n_studies)

# 2. processing and filtering raw output:
filtered_output <- raw_output |>
  process_output_nt() |>
  # retain only data for reliable gene sequences:
  filter_output_nt(min_seq_length = globsets$min_seq_length,
                    min_alig_score = globsets$min_alig_score,
                    max_core_dist = globsets$max_core_dist) |>
  # retain only data for mutations of interest:
  semi_join(mutation_list_reports, by = join_by(Nt_pos_Ecoli, Nt_mutation))

write_csv(filtered_output,"./output/rrs_filtered_output.csv")

# 3. analysis of extracted gene sequences and filtering:
summarise_target_sequences(genome_summaries, 
                           raw_output, 
                           filtered_output, 
                           min_seq_length = globsets$min_seq_length,
                           min_alig_score = globsets$min_alig_score,
                           max_core_dist = globsets$max_core_dist,
                           target_gene = 'rrs',
                           file_name = "./results/rrs_summary_target_sequences_rrs.txt")
plot_target_sequences_stats(raw_output, 
                            filtered_output,
                            min_seq_length = globsets$min_seq_length,
                            min_alig_score = globsets$min_alig_score,
                            max_core_dist = globsets$max_core_dist,
                            file_names = c("./plots/rrs_target_sequence_stats_hist.pdf", "./plots/rrs_target_sequence_stats_pairs.pdf"))

#empty working environment to keep everything clean:
rm.all.but("globsets")

########################################################################
### Step 6: Analysis of results  ###
########################################################################

# input for this step:  extracted target gene sequences (./"output/rrs_target_sequences.fa")
#                       E. coli gene reference sequence ("./data/rrs_references.fasta")
#                       filtered output of mutation screen ("./output/rrs_filtered_output.csv")
#                       bacterial taxonomic information ("./output/rrs_NCBI_taxonomy.csv")
#                       
# output for this step: filtered results for reliable sequences("./output/rrs_filtered_output.csv")
#                       summary of screening result ("./output/rrs_summary_screen_mutations.txt")
#                       plot of predictions for each mutation across species ("./plots/rrs_mutation_screening.pdf")
#                       plot of predictions for different classes and genera ("./plots/rrs_classes_genera.pdf")
#                       table of statistics for species with multiple gene copies ("./output/rrs_multiseq_stats.csv")
#                       plot showing the statistics for species with multiple gene copies ("./plots/rrs_multicopy_stats.pdf")

# 1. load required data:
rrs_target_sequences <- readDNAStringSet("./output/rrs_target_sequences.fa")
rrs_reference_Ecoli <- readDNAStringSet("./data/rrs_references.fasta")[["rrs_Escherichia_coli_MG1655"]]
filtered_output <- read_csv("./output/rrs_filtered_output.csv", show_col_types = FALSE)
bacterial_taxonomy <- read_csv("./data/rrs_NCBI_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from NCBI

#2. analysis of mutant screen:
plot_mutation_screen(filtered_output, file_name = "./plots/rrs_mutation_screen.pdf")
plot_classes_genera(filtered_output, bacterial_taxonomy, file_name= "./plots/rrs_classes_genera.pdf")
plot_evolvability_by_class(filtered_output, bacterial_taxonomy, file_name = "./plots/rrs_evolvability_by_class.pdf")
summarise_mutation_screen(filtered_output, target_gene = "rrs", file_name = "./results/rrs_summary_mutation_screen.txt")
get_resistance_taxonomy(filtered_output, bacterial_taxonomy, file_path = "./output/")
make_table_intrinsic_resistance(filtered_output, file_name = "./results/rrs_predicted_resistance.csv")

#3. analyse species with multiple gene copies:
multiseq_stats <- compare_gene_copies(filtered_output, rpsL_target_sequences, rpsL_reference_Ecoli)
write_csv(multiseq_stats, "./output/rrs_multiseq_stats.csv")
#multiseq_stats <- read_csv("./output/multiseq_stats.csv", show_col_types = FALSE)
plot_multiseq_stats(multiseq_stats, "./plots/rrs_multiseq.pdf")

#empty working environment to keep everything clean:
rm.all.but("globsets")


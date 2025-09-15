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
library(cowplot) 
library(ggraph)
library(tidygraph)
library(igraph)


source("R/util.R")
source("R/bioinformatics.R")
source("R/bioinformatics_nt.R")
source("R/analyses_nt.R")
source("R/analyses.R")
source("R/plotting_nt.R")
source("R/reports_nt.R")
source("R/reports.R")
source("R/phylogenetics.R")
source("R/structure.R")
source("R/process_nt_data.R")


# Global settings and parameters:
globsets <- list(
  min_n_studies = 2, # minimum number of studies that a mutation needs to be reported in for inclusion
  min_n_species = 2, # minimum number of species that a mutation needs to be reported in for inclusion
  min_seq_length = 1200, # minimum length of included gene target sequences
  min_alig_score = -Inf, # minimum alignment score (with E. coli) of included gene target sequences
  max_core_dist = 190, # maximum Levenshtein distance between E. coli core gene region to corresponding target region
  phylo_stats_sample_n = 5000, # number of species to sample for phylogenetics statistics
  random_seed = 22)
options(nwarnings = 10000)

set.seed(globsets$random_seed)

# Load data
meta_data <- read_tsv("data/rrnDB-5.10.tsv")
rrs_rrndb_sequences <- readDNAStringSet("./data/rrnDB-5.10_16S_rRNA.fasta")

# Define the desired database 
db <- "assembly" 
# define search term for representative genomes: Bacterial genomes at all assembly levels with annotation
term_rep <- '("Bacteria"[Organism] OR bacteria[All Fields]) AND ("latest refseq"[filter] AND "representative genome"[filter])'#representative
# search the entire database using the defined term
summaries_rep <- get_summaries(db, term_rep)
# define search term for reference genomes: Bacterial genomes at all assembly levels with annotation 
term_ref <- '("Bacteria"[Organism] OR bacteria[All Fields]) AND ("latest refseq"[filter] AND "reference genome"[filter])'
# search the entire database using the defined term
summaries_ref <- get_summaries(db, term_ref)
# combine summaries and save them
summaries <- c(summaries_ref, summaries_rep)

# Create assembly accession list
assembly_numbers_ref <- sapply(summaries, function(x) x$assemblyaccession) |> unname()

# Extract headers 
headers <- names(rrs_rrndb_sequences)
# Select sequences
selected_headers <- tibble(header = headers) |>
  separate(header,
           into = c("species", "assembly_number", "refseq_accession_number", "chromosome_info", "location"),
           sep = "\\|",
           remove = FALSE) |>
  filter(assembly_number %in% assembly_numbers_ref) |> 
  group_by(assembly_number) |>
  mutate(seq_index = row_number(),
         gene_copy = n()) |>
  ungroup() |>
  mutate(new_name = paste0("rrs_", seq_index, "_", species, "_", assembly_number))

# length(unique(selected_headers$assembly_number))
# 5896

# Subset sequences
selected_sequences <- rrs_rrndb_sequences[selected_headers$header]

# Apply new names
names(selected_sequences) <- selected_headers$new_name

# Save to fasta file
writeXStringSet(selected_sequences, "./output/rrs_target_sequences_rrnDB.fa")

#empty working environment to keep everything clean:
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
rrs_target_sequences <- readDNAStringSet("./output/rrs_target_sequences_rrnDB.fa")
rrs_reference_Ecoli <- readDNAStringSet("./data/rrs_reference_rrnDB.fasta")[["rrs_Escherichia_coli_U_5/41"]]

# 2.make a list of reliable mutations to be screened:
mutation_list_reports <- filter_mutations_nt(muts,
                                             min_n_species = globsets$min_n_species, 
                                             min_n_studies = globsets$min_n_studies)
mutation_list <- mutation_list_reports |>
  distinct() |>
  arrange(Nt_pos_Ecoli, Nt_mutation)

# 3.screen all rrs sequences for existing and possible mutations:
raw_output <- screen_target_sequences_nt(rrs_target_sequences, rrs_reference_Ecoli, 
                                         mutation_list, target_gene="rrs", n_workers=8)

#save error messages:
saveRDS(raw_output[!sapply(raw_output, is.data.frame)], "./output/rrs_raw_output_errors_rrnDB.rds")

#save results:
raw_output <- do.call(rbind, raw_output[sapply(raw_output, is.data.frame)])
write_csv(raw_output, file = "./output/rrs_raw_output_rrnDB.csv")

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
                           file_name = "./results/summary_rrs_target_sequences.txt")
plot_target_sequences_stats_nt(raw_output, 
                               filtered_output,
                               min_seq_length = globsets$min_seq_length,
                               min_alig_score = globsets$min_alig_score,
                               max_core_dist = globsets$max_core_dist,
                               file_names = c("./plots/rrs_target_sequence_stats_hist.pdf", "./plots/rrs_target_sequence_stats_pairs.pdf"))

#empty working environment to keep everything clean:
rm.all.but("globsets")













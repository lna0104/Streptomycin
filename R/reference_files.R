# Load necessary libraries
library(rentrez)
library(Biostrings)
library(tidyverse)

# Source additional R scripts for bioinformatics functions and utilities
source("R/bioinformatics.R")
source("R/util.R")

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
  
  file_path <- download_file(summaries[[j]], dir = "output/references")
  
  fasta_file <- sub("\\.gz$", "", file_path)
  sequences <- readDNAStringSet(fasta_file)
  matching_seq <- sequences[grep("rpsL|30S ribosomal protein S12($|\\])", names(sequences))]
  
  if (length(matching_seq) != 0) {
    names(matching_seq) <- added_name_muts[added_name_muts$ID == id, ]$FASTA_name
    all_sequences <- c(all_sequences, matching_seq)
  }else{
    cat(paste0("No rpsL found for species ID ", id, "\n"))
    next
  }
}
# Combine all sequences into a single DNAStringSet object
combined_sequences <- do.call(c, all_sequences)

# Write the combined sequences to a new FASTA file
writeXStringSet(combined_sequences, "data/rpsL_references.fasta")





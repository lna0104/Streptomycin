
#' Download summaries from an Entrez database search
#'
#' Searches an Entrez database for a given search term, downloading the summaries for each result.
#'
#' @param db a character string giving the name of the Entrez database to search.
#' @param term a character string giving the search term to pass to the database.
#' @param retmax the maximum number of summaries to request from the database at a time.
#' @param max_tries the maximum number of times to try each web request before throwing an error.
#'
#' @details
#' `retmax` limits the number of summaries requested in each web request to the database.
#' `get_summaries()` will make multiple requests, each of size not exceeding `retmax`, until all available summaries have been downloaded.
#' Each web request will be attempted `max_tries` times to protect against intermittent internet failure.
#'
#' @returns a list of esummary records
#' @export
#'
#' @examples
#' get_summaries('assembly', 'GCA_000006945.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])')
get_summaries <- function(db, term, retmax = 500, max_tries = 10) {
  #Firstly, we perform a search, but set retmax = 0 to have it return no entries.
  #This just loads entries into the history server.
  #We will then use web_history to get the summaries.
  search_result <- rentrez::entrez_search(db = db, term = term, retmax = 0, use_history = TRUE)
  #Entrez limits the number of entries it will return in a single search.
  #As such, we need to make repeated requests to get all the summaries.
  summaries <- list()
  while(length(summaries) < search_result$count) {
    #The web request may fail, and we don't want to lose all our progress so far if it does, so use retry() for it.
    these_summaries <- retry(
      rentrez::entrez_summary(db = db, web_history = search_result$web_history, retstart = length(summaries), retmax = retmax, always_return_list = TRUE),
      max_tries = max_tries
    )
    summaries <- append(summaries, these_summaries)
  }
  return(summaries)
}

# extracts the file path to a genome from a summary:
get_genome_file_path <- function(summary, suffix, dir) {
  name <- gsub(".*/([^/]+)/?", "\\1", summary$ftppath_refseq)
  file_name <- paste0(name, suffix)
  file_path <- file.path(dir, file_name)
  return(file_path)
}


#Takes an Entrez summary, and downloads and unzips the fasta file from refseq.
#Each summary object must have an `ftppath_refseq` entry.
#Downloaded files are placed in the `dir` directory, which is created if it does not exist.
#' Download Fasta File From RefSeq
#'
#' Given an esummary object, downloads (and unzips if necessary) a corresponding fasta files of a given type from RefSeq.
#'
#' @param summary an esummary object with an `ftppath_refseq` entry.
#' @param suffix a character string giving the suffix of the file to be downloaded (following the assembly name).
#' @param dir a character string giving the directory to place the downloaded file in (created if it does not exist).
#' @param max_tries the maximum number of times to try each web request before throwing an error.
#'
#' @details
#' The URL from which to download the file is constructed from `summary`'s `ftppath_refseq` entry and the given `suffix`.
#' `suffix` must include the file extension, including `.gz` if applicable.
#' Each web request will be attempted `max_tries` times to protect against intermittent internet failure.
#' The downloaded file is unzipped if it ends in `gz`, and the original `.gz` files is removed.
#'
#' @returns a character string giving the path to the downloaded (and unzipped if applicable) file
#' @export
#'
#' @examples
#' withr::with_tempfile("td", {
#'   summaries <- get_summaries('assembly', 'GCA_000006945.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])')
#'   files <- list()
#'   for(i in 1:length(summaries)) {
#'     try( #Continue if any download fails
#'       files[[i]] <- download_file(summaries[[i]], dir = td)
#'     )
#'   }
#' })
download_file <- function(summary, suffix = "_cds_from_genomic.fna.gz", dir = "output/genomes", max_tries = 10, overwrite = FALSE) {
  #Create the output directory if it does not exist.
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  name <- gsub(".*/([^/]+)/?", "\\1", summary$ftppath_refseq)
  file_name <- paste0(name, suffix)
  # file_path <- file.path(dir, file_name)
  file_path <- get_genome_file_path(summary, suffix = suffix, dir = dir)

  if ((!overwrite) && (file.exists(sub(".gz", "", file_path)))) {
    cat("  --> File already exists, skipped.\n")
  } else {
     url <- paste0(summary$ftppath_refseq, "/", file_name)
     url <- gsub("^ftp:", "https:", url)
     print(url)
     retry(
       download.file(url, file_path),
       max_tries = max_tries
     )
     if(grepl("gz$", suffix)) {
       file_path <- R.utils::gunzip(file_path, overwrite = TRUE)
     }
  }
  return(file_path)
}

#Takes an list of Entrez summaries, and downloads and unzips the fasta files from refseq.
#Each summary object must have an `ftppath_refseq` entry.
#Downloaded files are placed in the `dir` directory, which is created if it does not exist.
#' Download Fasta Files From RefSeq
#'
#' Given a list of esummary objects, downloads (and unzips if necessary) a set of corresponding fasta files of a given type from RefSeq.
#'
#' @param summaries a list of esummary objects with a `ftppath_refseq` entries.
#' @param suffix a character string giving the suffix of the files to be downloaded (following the assembly name).
#' @param dir a character string giving the directory to place the downloaded files in (created if it does not exist).
#' @param max_tries the maximum number of times to try each web request before throwing an error.
#'
#' @details
#' The URLs from which to download the files are constructed from the `ftppath_refseq` entry of each summary and the given `suffix`.
#' `suffix` must include the file extension, including `.gz` if applicable.
#' Each web request will be attempted `max_tries` times to protect against intermittent internet failure.
#' The downloaded files are unzipped if they end in `gz`, and the original `.gz` files are removed.
#'
#' @returns a vector of character strings giving the paths to the downloaded (and unzipped if applicable) files
#' @export
#'
#' @examples
#' withr::with_tempfile("td", {
#'   summaries <- get_summaries('assembly', 'GCA_000006945.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])')
#'   files <- download_files(summaries, dir = td)
#' })
download_files <- function(summaries, suffix = "_cds_from_genomic.fna.gz", dir = "output/genomes", max_tries = 10) {
  file_paths <- list()
  for(i in 1:length(summaries)) {
    cat(paste0("Downloading genome for species #", i, "/", length(summaries), "\n..."))
    try( # continue to download the genomes if any download fails
      file_path <- download_file(summaries[[i]], suffix = suffix, dir = dir, max_tries = max_tries)
    )
    file_paths[[i]] <- file_path
  }
  return(file_paths)
}



#' Find the target sequence inside the annotated genomic file based on the presence of the gene or protein name 
#'
#' @param fna_file a genomic(cds) fna file
#' @param target_gene a character string object giving the name of a gene("rpoB") to be retrieved
#' Find all rpoB sequences inside the annotated genomic file based on the presence of the gene or protein name 
#'
#' @param fna_file a genomic(cds) fna file
#' @param target_gene a character string object giving the name of a gene("rpoB") to be retrieved
#' @param target_protein a character string object giving the name of a protein ("rpoB") to be retrieved. 
#' This is used as an alternative to target_gene.
#'
#' @return DNAStringSet object containing the genes of interest 
#' @export
#'
#' @examples get_target_genes(cds_from_genome.fna)
get_target_genes <- function(fna_file, target_gene, target_protein) {
  lines <- readLines(fna_file)
  
  # We want to split the lines into a list of headings
  all_genes <- base::split(lines, cumsum(grepl("^>", lines)))
  target_gene_list <- list()
  for(gene in all_genes) {
    if (grepl(sprintf("\\[gene=%s\\]", target_gene), gene[1]) || 
        grepl(sprintf("\\[protein=%s\\]", target_protein), gene[1])) {
      # Stitch all the lines together, except the first
      target_gene_list <- c(target_gene_list, DNAString(paste(gene[-1], collapse = "")))
    }
  }
  return(target_gene_list)
}


#' Extract the target sequences from a set of downloaded genomes
#'
#' @param summaries A list of multiple esummaries from an entrez database containing information about genomes
#' @param file_paths A list of paths to downloaded fasta files, corresponding to the list of summaries (not necessarily in the same order)
#' @param target_gene a character string object giving the name of a gene ("rpoB") to be retrieved
#' @param target_protein a character string object giving the name of a protein ("rpoB") to be retrieved. 
#'
#' @return a list of DNA string objects (rpoB sequences) and their corresponding identifications (species name and accession number), in the same order as the summaries
#' @export 
#'
#' @examples
#' withr::with_tempfile("td", {
#'   summaries <- get_summaries('assembly', 'GCA_000006945.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])')
#'   files <- download_files(summaries, dir = td)
#'   target_sequences <- get_target_sequences(summaries, files)
#' })
get_target_sequences <- function(summaries,
                                      suffix = "_cds_from_genomic.fna", 
                                      dir = "output/genomes",
                                      target_gene = 'rpoB',
                                      target_protein = "DNA-directed RNA polymerase subunit beta($|[^'])|DNA-directed RNA polymerase subunit beta chain|RNA polymerase, beta subunit|DNA-directed RNA polymerase subunit beta/beta'") {
  # file_paths <- sort(unlist(file_paths))
  # file_paths <- file_paths[order(order(unlist(lapply(summaries, `[[`, "assemblyaccession"))))]
  target_sequences <- list()
  for(i in 1:length(summaries)) {
    cat(paste0("Extracting gene(s) for species #", i, "/", length(summaries), "\n..."))
    file_path <- get_genome_file_path(summaries[[i]], suffix = suffix, dir = dir)
    targets <- get_target_genes(file_path, 
                                target_gene = target_gene, 
                                target_protein = target_protein) #extract the rpoB genes
    if (length(targets) == 0L) {
      targets <- list(DNAString())  
    }
    names(targets) <- paste0(target_gene,"_", 1:length(targets), "_", summaries[[i]]$speciesname, "_", summaries[[i]]$assemblyaccession)
    target_sequences <- c(target_sequences, targets)
  }
  return(target_sequences)
}



#' Download taxonomy data from NCBI
#'
#' @param summaries A list of multiple esummaries from an entrez database containing information about genomes
#' @param database_file Name of the file for the taxonomy database containing all taxonomy data from NCBI
#' @param output_file Name of the file to be created containing taxonomy data for species in summary
#'
#' @return NULL
#' @export
#'
download_taxonomy <- function(summaries,
                              temp_folder = "./data/tmp/",
                              output_file = "./data/NCBI_taxonomy.csv") {
  
  if (!dir.exists(temp_folder)) {
    dir.create(temp_folder)
  }
  database_file = paste0(temp_folder, "nameNode.sqlite") 
  # download the database:
  taxonomizr::prepareDatabase(database_file, tmpDir = "./data/taxonomizr-database", getAccessions = FALSE, protocol = "http")
  
  our_taxa <- data.frame(tax_id = rep(NA, length(summaries)),
                         species = rep(NA, length(summaries)))
  for (i in 1:length(summaries)){
    our_taxa$tax_id[i] <- summaries[[i]]$taxid
    our_taxa$species[i] <- summaries[[i]]$speciesname
  }
  our_taxa <- cbind(our_taxa, 
                    taxonomizr::getTaxonomy(our_taxa$tax_id, 
                                            sqlFile = database_file, 
                                            desiredTaxa = c("phylum", "class", "order", "family", "genus"))) |>
    select(-tax_id, -species) |>
    distinct() |>
    # add class to two genera where class is missing:
    mutate(class = ifelse(genus == "Mycoplasmopsis", "Mollicutes", class)) |>
    mutate(class = ifelse(genus == "Caldicellulosiruptor", "Clostridia", class))

  write_csv(our_taxa, output_file)
  # delete temporary files:
  unlink(temp_folder, recursive = TRUE)
  return(invisible(NULL))
}

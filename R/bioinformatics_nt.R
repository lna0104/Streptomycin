#' Extract the target sequences from a set of downloaded genomes
#'
#' @param summaries A list of multiple esummaries from an entrez database containing information about genomes
#' @param file_paths A list of paths to downloaded fasta files, corresponding to the list of summaries (not necessarily in the same order)
#' @param target_gene a character string object giving the name of a gene ("rrs") to be retrieved
#' @param target_rna a character string object giving the name of a RNA to be retrieved. 
#'
#' @return a list of DNA string objects (rrs sequences) and their corresponding identifications (species name and accession number), in the same order as the summaries
#' @export 
#'
#' @examples
#' withr::with_tempfile("td", {
#'   summaries <- get_summaries('assembly', 'GCA_000006945.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])')
#'   files <- download_files(summaries, dir = td)
#'   target_sequences <- get_target_sequences(summaries, files)
#' })
get_target_sequences_nt <- function(summaries,
                                      suffix = "_cds_from_genomic.fna", 
                                      dir = "output/genomes",
                                      target_gene,
                                      target_rna) {
  # file_paths <- sort(unlist(file_paths))
  # file_paths <- file_paths[order(order(unlist(lapply(summaries, `[[`, "assemblyaccession"))))]
  target_sequences <- list()
  for(i in 1:length(summaries)) {
    cat(paste0("Extracting gene(s) for species #", i, "/", length(summaries), "\n..."))
    file_path <- get_genome_file_path(summaries[[i]], suffix = suffix, dir = dir)
    targets <- get_target_genes_nt(file_path, 
                                target_gene = target_gene, 
                                target_rna = target_rna) #extract the target genes
    if (length(targets) == 0L) {
      targets <- list(DNAString())  
    }
    names(targets) <- paste0(target_gene,"_", 1:length(targets), "_", summaries[[i]]$speciesname, "_", summaries[[i]]$assemblyaccession)
    target_sequences <- c(target_sequences, targets)
  }
  return(target_sequences)
}

#' Find the target sequence inside the annotated genomic file based on the presence of the gene or rna name 
#'
#' @param fna_file a genomic(cds) fna file
#' @param target_gene a character string object giving the name of a gene("rpoB") to be retrieved
#' Find all rpoB sequences inside the annotated genomic file based on the presence of the gene or protein name 
#'
#' @param fna_file a genomic(cds) fna file
#' @param target_gene a character string object giving the name of a gene("rpoB") to be retrieved
#' @param target_rna a character string object giving the name of a RNA to be retrieved. 
#' This is used as an alternative to target_gene.
#'
#' @return DNAStringSet object containing the genes of interest 
#' @export
#'
#' @examples get_target_genes(rna_from_genome.fna)
get_target_genes_nt <- function(fna_file, target_gene, target_rna) {
  lines <- readLines(fna_file)
  
  # We want to split the lines into a list of headings
  all_genes <- base::split(lines, cumsum(grepl("^>", lines)))
  target_gene_list <- list()
  for(gene in all_genes) {
    if (grepl(sprintf("\\[gene=%s\\]", target_gene), gene[1]) || 
        grepl(sprintf("\\[product=%s\\]", target_rna), gene[1])) {
      # Stitch all the lines together, except the first
      target_gene_list <- c(target_gene_list, DNAString(paste(gene[-1], collapse = "")))
    }
  }
  return(target_gene_list)
}

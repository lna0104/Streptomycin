# List of packages you loaded
pkgs <- c(
  "ALJEbinf", "tidyverse", "openssl", "stringr", "parallel", "future.apply", "varhandle", "rentrez", "Biostrings", "pwalign", "MSA2dist", "castor", "ape", "phytools", "phangorn", "ggtree", "tidytree", "treeio", "bio3d", "ggnewscale", "colorspace", "scales", "patchwork", "ggpubr", "plotrix", "RColorBrewer", "GGally", "ggh4x",
  "pander", "quarto", "NGLVieweR", "htmlwidgets", "wesanderson", "ggrepel", "cowplot", "ggraph", "tidygraph", "igraph"
)

# Function to get BibTeX citations for each package
get_bibtex_citations <- function(pkgs) {
  bibs <- lapply(pkgs, function(pkg) {
    # Try to get citation, some packages might not have citation info or installed
    tryCatch({
      # citation() might return multiple citations, convert all
      cite <- citation(pkg)
      bibtex_entries <- sapply(cite, toBibtex)
      return(bibtex_entries)
    }, error = function(e) {
      message(sprintf("No citation found for package: %s", pkg))
      return(NULL)
    })
  })
  # Flatten and remove NULLs
  bibs <- unlist(bibs)
  bibs <- bibs[!sapply(bibs, is.null)]
  return(bibs)
}

# Get BibTeX citations
bibtex_citations <- get_bibtex_citations(pkgs)

# Write to file
output_file <- "./output/packages_citations.bib"
writeLines(bibtex_citations, con = output_file)

cat(sprintf("BibTeX citations saved to '%s'\n", output_file))

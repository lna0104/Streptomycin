library(Biostrings)
library(msa)

rrs_reference_Ecoli_new <- readDNAStringSet("./data/rrs_reference_rrnDB.fasta")[["rrs_Escherichia_coli_U_5/41"]]
rrs_reference_Ecoli_old <- readDNAStringSet("./data/rrs_references.fasta")[["rrs_Escherichia_coli_MG1655"]]

seqs <- DNAStringSet(c(
  old = rrs_reference_Ecoli_old,
  new = rrs_reference_Ecoli_new
))

aln <- msa(seqs)

msaPrettyPrint(
  aln,
  output = "pdf",
  showNames = "left",     # show sequence names
  showLogo = "none",      # no sequence logo (optional)
  askForOverwrite = FALSE,
  file = "plots/alignment_msa.pdf"  # output filename
)

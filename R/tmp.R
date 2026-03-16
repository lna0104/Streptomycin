library(dplyr)
library(tidyr)
library(stringr)
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

library(ggplot2)
library(ggpubr)

source("R/analyses.R")
source("R/util.R")


# Functions
filter_mutations <- function(muts,
                             min_n_species = 1,
                             min_n_studies = 1,
                             origin = "ANY") {
  mutation_list <- muts |>
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_original), !is.na(AA_mutation))
  # filter(is.na(Warning) | Warning == "AA_pos inconsistent with AA_pos_Ecoli")

  if (origin == "ANY") {
    origin <- muts |>
      pull(Origin) |>
      unique()
  }

  mutation_list_species <- mutation_list |>
    select(Species, AA_pos_Ecoli, AA_mutation) |>
    distinct() |>
    group_by(AA_pos_Ecoli, AA_mutation) |>
    summarise(n_species = n(), .groups = "drop")

  mutation_list_studies <- mutation_list |>
    select(Ref_code, AA_pos_Ecoli, AA_mutation) |>
    distinct() |>
    group_by(AA_pos_Ecoli, AA_mutation) |>
    summarise(n_studies = n(), .groups = "drop")

  mutation_list <- mutation_list_studies |>
    full_join(mutation_list_species, by = join_by(AA_pos_Ecoli, AA_mutation)) |>
    filter(n_studies >= min_n_studies) |>
    filter(n_species >= min_n_species) |>
    select(AA_pos_Ecoli, AA_mutation)

  return(mutation_list)
}

check_out_mutation <- function(mutations_list,
                               target_sequence,
                               reference_Ecoli,
                               target_gene,
                               alig_file = NULL) {
  if (is.null(target_sequence) || length(target_sequence) == 0) {
    stop(sprintf('Empty "%s" sequence.', target_gene))
  }
  coords_alig <- ALJEbinf::getCoordinates(target_sequence, reference_Ecoli, aligOutput = TRUE)
  if (!is.null(alig_file)) {
    save(coords_alig, file = alig_file)
  }
  coords <- coords_alig$coordinates
  protein_target <- Biostrings::translate(target_sequence)
  protein_Ecoli <- Biostrings::translate(reference_Ecoli)

  result <- mutations_list
  result$AA_target <- NA # what is the amino acid at that position in focal species
  result$Codon_target <- NA # what's the codon at that position in focal species
  result$n_possible <- NA # number of single nt mutations in sequence that would produce the AA mutation
  for (i in 1:nrow(mutations_list)) {
    AA_pos <- ALJEbinf::translateCoordinate(mutations_list$AA_pos_Ecoli[i],
      coords,
      direction = "RefToFocal",
      AAinput = TRUE,
      AAoutput = TRUE
    )
    AA_mutation <- mutations_list$AA_mutation[i]
    if (!is.na(AA_pos)) { # only proceed if corresponding AA presents in target sequence!
      result$AA_target[i] <- as.character(protein_target[AA_pos])
      nt_pos <- AA_pos * 3 - 2
      result$Codon_target[i] <- as.character(target_sequence[nt_pos:(nt_pos + 2)])
      result$n_possible[i] <- n_mutations_in_codon(
        AA_pos,
        mutations_list$AA_mutation[i],
        target_sequence
      )
    }
  }
  result$alig_score <- pwalign::score(coords_alig$alignment) # alignment score for each target sequence

  # calculate Levenshtein distance between E. coli and target gene core region:
  # determine whether to use the core region or the original sequence

  from_target <- ALJEbinf::translateCoordinate(40, coords, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
  to_target <- ALJEbinf::translateCoordinate(100, coords, direction = "RefToFocal", AAinput = TRUE, AAoutput = TRUE)
  cat(paste0("Working on ", from_target, "_", to_target, "...\n"))
  core_target <- protein_target[from_target:to_target]
  core_Ecoli <- protein_Ecoli[40:100]

  result$core_dist <- pwalign::stringDist(c(Biostrings::AAStringSet(core_target), Biostrings::AAStringSet(core_Ecoli)),
    method = "levenshtein"
  )
  return(result)
}

screen_target_sequences <- function(target_sequences, reference_Ecoli, mutation_list, target_gene, n_workers) {
  # Ensure required directories exist
  if (!dir.exists("./output/alignments_rickettsia")) dir.create("./output/alignments_rickettsia", recursive = TRUE)
  if (!dir.exists("./output/progress_rickettsia")) dir.create("./output/progress_rickettsia", recursive = TRUE)

  # plan(multisession, workers = n_workers) # to parallelise the function
  final_output <- future_lapply(1:length(target_sequences), function(i) {
    cat(paste0("Working on species #", i, "/", length(target_sequences), "...\n")) # notify when the analysis of each sequence starts
    output <- tryCatch(
      check_out_mutation(mutation_list, # compare target sequences with reference sequence to identify the amino acid at each position
        target_sequences[[i]],
        reference_Ecoli,
        target_gene,
        alig_file = paste0("./output/alignments_rickettsia/", names(target_sequences)[[i]], ".RData")
      ),
      error = function(cond) {
        message("  It seems there is a problem with this sequence. Here's the original error message:") # gives error but keeps the analysis moves on
        message(conditionMessage(cond)) # show what is the error
        return(list(error = conditionMessage(cond)))
      }
    )
    if (!is.null(output)) {
      output$accession_numbers <- extract_accession_number(names(target_sequences)[i], target_gene) # insert a column for each target gene assembly accession number
      output$species <- extract_species_name(names(target_sequences)[i], target_gene) # insert a column for each species name
      output$target_name <- names(target_sequences)[i]
      output$gene_copy <- extract_gene_copy_number(names(target_sequences)[i], target_gene) # insert a column for each target gene copy number
      output$target_length <- length(target_sequences[[i]])
      if (is.data.frame(output)) {
        output <- select(output, accession_numbers:target_length, alig_score, core_dist, everything())
      }
    }

    saveRDS(i, file = paste0("./output/progress_rickettsia/sequence_", i, "_out_of_", length(target_sequences), "_complete.rds"))
    return(output)
  })
  return(final_output)
}

# Analysis

muts <- read.csv("./output/rpsL_checked_muts.csv")
raw_output <- read_csv("./output/rpsL_raw_output.csv", show_col_types = FALSE)
rpsL_target_sequences <- readDNAStringSet("./output/rpsL_target_sequences.fa")
rpsL_reference_Ecoli <- readDNAStringSet("./data/rpsL_references.fasta")[["rpsL_Escherichia_coli_MG1655"]]
rickettsia_hits <- names(rpsL_target_sequences)[grepl("Rickettsia", names(rpsL_target_sequences))]
rickettsia_seqs <- rpsL_target_sequences[grepl("Rickettsia", names(rpsL_target_sequences))]


# 2.make a list of reliable mutations to be screened:
mutation_list_reports <- filter_mutations(muts,
  min_n_species = 3,
  min_n_studies = 3
)

mutation_list <- mutation_list_reports |>
  distinct() |>
  arrange(AA_pos_Ecoli, AA_mutation)

raw_output <- screen_target_sequences(rickettsia_seqs, rpsL_reference_Ecoli,
  mutation_list,
  target_gene = "rpsL", n_workers = 1
)

raw_output <- do.call(rbind, raw_output[sapply(raw_output, is.data.frame)])
write_csv(raw_output, file = "./output/rpsL_raw_output_rickettsia.csv")

filtered_output_rickettsia <- raw_output |>
  process_output() |>
  semi_join(mutation_list, by = join_by(AA_pos_Ecoli, AA_mutation))

load("./output/alignments_rickettsia/rpsL_1_Rickettsia helvetica_GCF_963970025.1.RData")
View(coords_alig[["coordinates"]])


dat_raw <- select(raw_output, species, accession_numbers, gene_copy, target_length, alig_score, core_dist) |>
  mutate(filter = factor("raw")) |>
  distinct()

# Plotting
plot1 <- ggplot(dat_raw) +
  geom_histogram(aes(target_length), position = "identity") +
  labs(x = "Length", y = "Number of sequences", fill = "")
plot2 <- ggplot(dat_raw) +
  geom_histogram(aes(alig_score), position = "identity", binwidth = 20) +
  labs(x = "Alignment score", y = NULL, fill = "")
plot3 <- ggplot(dat_raw) +
  geom_histogram(aes(core_dist), position = "identity", binwidth = 3) +
  labs(x = "Levenshtein distance", y = NULL, fill = "")
hists <- ggarrange(plot1, plot2, plot3,
  nrow = 1,
  labels = LETTERS[1:3],
  common.legend = TRUE, vjust = 0
)
ggsave(filename = "./plots/rpsL_target_sequence_stats_hist_rickettsia.pdf", hists, width = 10, height = 4)

# rrs Literature Review
muts <- read.csv("./data/reported_mutations.csv") |>
  filter(Gene == "rrs") |>
  filter(Ref_code != "Viet2018")

ref_muts <- muts |>
  distinct(Species, Origin, Ref_code, Evidence_category)

p <- ggplot(ref_muts, aes(x = Evidence_category)) +
  geom_bar() +
  theme_minimal() +
  labs(
    x = "Evidence category",
    y = "Number of papers"
  )

ggsave(
  "./plots/rrs_evidence_category_literature.pdf",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

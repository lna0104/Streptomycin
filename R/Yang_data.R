# Functions to import mutations from Yang et al. (2023) paper


#' Import and wrangle data from Yang et al. (2023)
#'
#' @param file_name_RIF name of file containing frequencies after RIF exposure (Fig.S1b)
#' @param file_name_t0 name of file containing frequencies before RIF exposure (Fig. S1a)
#' @param min_log2_fc minimum log2(fold-change) of frequencies for mutations to be included 
#' as RIF-resistance mutations
#'
#' @return A table of mutations from Yang et al. that are classified as RIF-resistance mutations
#'
#' @examples
import_Yang_muts <- function(file_name_RIF_log2FC, min_log2_fc) {
  mut_list_Yang <- read_csv(file_name_RIF_log2FC, show_col_types = FALSE) |>
    pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log2_fc") |>
    filter(log2_fc >= min_log2_fc) |>
    mutate(pos = as.integer(paste0("5", pos))) |>
    rename(AA_pos_Ecoli = pos) |>
    select(AA_pos_Ecoli, AA_mutation)
  return(mut_list_Yang)
}

#' Plot comparison between mutations reported in the literature vs. in Yang et al. (2023)
#'
#' @param mut_list_reported Table of mutations reported as RIF-resistance mutations in the literature.
#' @param mut_list_Yang Table of mutations reported as RIF-resistance mutations in Yang et al. (2023)
#' @param file_name Name of pdf file for the plot
#'
#' @return NULL. Saves a pdf file with the plot.
#'
#' @examples
plot_reported_vs_Yang <- function(mutation_list_reports, mutation_list_Yang, file_name) {
  all_AAs_sorted <- str_split_1("HAKRGILPVFWYDESTCMNQ", "")
  mutation_comparison <- mutation_list_reports |>
    mutate(Source = "reported") |>
    rbind(mutation_list_Yang |> mutate(Source = "Yang")) |>
    group_by(AA_pos_Ecoli, AA_mutation) |>
    reframe(n = n(),
            Comb_source = case_when(n() == 1 ~ Source,
                                    n() == 2 ~ "both")) |>
    mutate(AA_mutation = factor(AA_mutation, levels = all_AAs_sorted)) |>
    complete(AA_pos_Ecoli, AA_mutation) |>
    mutate(Comb_source = factor(ifelse(is.na(Comb_source), "none", Comb_source),
                                levels = c("reported", "Yang", "both", "none"))) |>
    mutate(AA_pos_Ecoli = factor(AA_pos_Ecoli))
  
  p <- ggplot(mutation_comparison) +
    geom_tile(aes(x = AA_mutation, y = AA_pos_Ecoli, fill = Comb_source)) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values = c("blue", "red", "purple", "grey90")) +
    labs(x = "Amino acid", y = "Position", fill = "") +
    theme_bw()
  
  ggsave(file_name, p, width = 5, height = 6)
  return(invisible(NULL))
}



check_Yang_muts <- function(file_name_RIF, file_name_t0, file_name_RIF_log2FC) {
  Yang_muts_RIF <- read_csv(file_name_RIF, show_col_types = FALSE) |>
    pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log10_rel_freq_RIF")
  Yang_muts_t0 <- read_csv(file_name_t0, show_col_types = FALSE) |>
    pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log10_rel_freq_t0")
  Yang_muts_RIF_log2FC_reported <- read_csv(file_name_RIF_log2FC, show_col_types = FALSE) |>
    pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log2_FC_reported")
  Yang_muts <- Yang_muts_RIF |>
    left_join(Yang_muts_t0) |>
    left_join(Yang_muts_RIF_log2FC_reported) |>
    mutate(log2_FC_calculated = log2((10^log10_rel_freq_RIF) / (10^log10_rel_freq_t0)))
  
  pow10 <- function(x) { 10^x }
  cat("Total frequency of relative mutant abundances with RIF and at t=0:\n")
  Yang_muts |> pull(log10_rel_freq_RIF) |> pow10() |> sum(na.rm = TRUE) |> print()
  Yang_muts |> pull(log10_rel_freq_t0) |> pow10() |> sum(na.rm = TRUE) |> print()
  
  print(plot(Yang_muts$log2_FC_calculated, Yang_muts$log2_FC_reported))
  abline(h=0, col = "red")
  abline(v=0, col = "red")
  cat("Linear model for reported vs. calculated log2(FC) values:\n")
  summary(lm(Yang_muts$log2_FC_reported ~ Yang_muts$log2_FC_calculated))
  
  # determine mutations with FC>1:
  selected_Yang_muts_reported <- Yang_muts |>
    filter(log2_FC_reported > 0) |>
    select(pos, AA_mutation)
  selected_Yang_muts_calculated <- Yang_muts |>
    filter(log2_FC_calculated > 0) |>
    select(pos, AA_mutation)
  cat("Mutations with calculated FC>1 but not reported FC>1:\n")
  anti_join(selected_Yang_muts_calculated, selected_Yang_muts_reported) |> print()
  cat("Mutations with reported FC>1 but not calculated FC>1:\n")
  anti_join(selected_Yang_muts_reported, selected_Yang_muts_calculated) |> print()
}



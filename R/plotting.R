# This file contains functions for producing plots and figures.


#' Produces figure giving the frequency of reported position, mutation and species based on Google sheet
#' @param muts a data frame table providing info on all reported mutations in the literature
#' @param file_name the path that the plot should be save in
#' @param n_frequency number of reported mutations, positions and species
#' @return Three plots as ggplot objects in one single pdf file
#' @export
#'
#' @examples plot_reported_mutations(muts, "./plot/myFilename.pdf")
plot_reported_mutations <- function(muts, file_name, n_frequency) {

  #1. the frequency of reported mutations across sites (currently not included in figure)
  frequency_position <-  muts %>% 
    select(Species, AA_pos_Ecoli, Origin) %>%
    filter(!is.na(AA_pos_Ecoli)) %>%
#    filter(AA_pos_Ecoli %in% mutation_list_reports$AA_pos_Ecoli) %>% 
    group_by(Species, AA_pos_Ecoli) %>%
    summarise(Origin = ifelse(all(Origin == "Isolate"), "Isolate",
                              ifelse(all(Origin == "Lab-generated"), "Lab-generated", "Both")),
              .groups = "drop")

  plot1 <- frequency_position %>%
    mutate(AA_pos_Ecoli = factor(AA_pos_Ecoli, levels = sort(unique(AA_pos_Ecoli)))) %>%
    ggplot() +
    geom_bar(aes(x = AA_pos_Ecoli, fill = Origin), width = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.1, 
                                     size = 10),
          axis.text.y = element_text(size = 10), 
          axis.title=element_text(size=10, face="bold"),
          plot.title=element_text(hjust=0),
          legend.position = "top")+
    labs(x = "Amino acid position",
         y = "Number of species", 
         fill = "Origin:") +
    scale_y_continuous(breaks = seq(0, 50, by = 5), expand = c(0.01, 0)) +
    scale_fill_brewer(palette = "Pastel1") 

  #2. the frequency of each reported amino acid substitution across species:

  frequency_mutation <- muts |>
    select(Species, AA_pos_Ecoli, AA_mutation, Origin) |>
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) |>
    mutate(mutation_name = paste0(AA_pos_Ecoli, AA_mutation)) |> 
    select(Species, mutation_name, Origin) %>%
    group_by(Species, mutation_name) %>%
    summarise(Origin = ifelse(all(Origin == "Isolate"), "Isolate",
                              ifelse(all(Origin == "Lab-generated"), "Lab-generated", "Both")),
              .groups = "drop") |>
    arrange(desc(mutation_name)) |>
    mutate(mutation_name = factor(mutation_name, levels = sort(unique(mutation_name)))) |>
    group_by(mutation_name) |>
    filter(n() >= n_frequency)
  
  plot2 <- ggplot(frequency_mutation) +
    geom_bar(aes(x = mutation_name, fill = Origin))  +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.1, 
                                     hjust=0.95, 
                                     size=6), 
          axis.text.y = element_text(size=6), 
          axis.title=element_text(size=6, face="bold"),
          plot.title=element_text(hjust=0),
          legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.position = "top") +
    labs(x = "Amino acid substitution",
         y = "Number of species", 
         fill = "Origin:") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_fill_brewer(palette = "Pastel1", breaks=c('Lab-generated', 'Isolate', 'Both')) 

  #3. how many mutations for each species have been reported
  frequency_per_species <-  muts |>
    select(Species, AA_pos_Ecoli, AA_mutation, Origin) |>
    filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation)) |>
    mutate(mutation_name = paste0(AA_pos_Ecoli, AA_mutation)) |>
    group_by(mutation_name, Species) |>
    summarise(Origin = ifelse(all(Origin == "Isolate"), "Isolate",
                              ifelse(all(Origin == "Lab-generated"), "Lab-generated", "Both")),
              .groups = "drop") |>
    group_by(Species) |> 
    filter(n() > n_frequency) 
  
  plot3 <- frequency_per_species %>%
    ggplot() +
    geom_bar(aes(x = fct_infreq(Species), fill = Origin), width = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.95,
                                     vjust=0.2, 
                                     size=6), 
          axis.title=element_text(size=6, face="bold"),
          axis.text.y=element_text(size=6), 
          plot.title=element_text(hjust=0), 
          legend.position = "top") +
    labs(x = "Species",
         y = "Number of reported mutations", 
         fill = "Origin:") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_fill_brewer(palette = "Pastel1") 

  #combine plots:
  reported_mutations <- ggarrange(plot2, plot3, 
                                  nrow=2, common.legend = TRUE,
                                  labels = c("A", "B"), hjust = 0, vjust = 0.5)
  
  ggsave(filename = file_name, reported_mutations, width = 6, height = 6)
}
  

#' Produces figure giving the frequency of mutations in three present, possible and impossible categories
#' @param filtered_output  a data frame providing screened and filtered gene sequences
#' @param file_name the path that the plot should be save in
#' @return A plot as a ggplot object in pdf format
#' @export
#'
#' @examples mutation_screen(filtered_output, "myFilename.pdf")
plot_mutation_screen <- function(filtered_output, file_name) {
  
  # species with more than one gene copy, to be filtered out:
  filtered_output_seqs <- filtered_output |>
    select(accession_numbers, species, gene_copy, target_length, alig_score, core_dist) |>
    distinct()
  multicopy_species <- filtered_output_seqs |>
    group_by(accession_numbers) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n > 1L) |>
    pull(accession_numbers)
  
  data_for_plotting <- filtered_output |>
    filter(!(accession_numbers %in% multicopy_species)) |>
    mutate(mutation_category = ifelse(mutation_category == "possible", 
                                      n_possible, 
                                      mutation_category)) |>
    mutate(mutation_category = factor(mutation_category, 
                                      levels = c("impossible", "1", "2", "3", "present"))) |>
    mutate(mutation_name = str_replace(mutation_name, "_", "")) |>
    mutate(mutation_name = factor(mutation_name, levels = sort(unique(mutation_name))))
    
  cols = c(hsv(0, 0, c(0.2, 0.73, 0.8)), hsv(0, 1, 0.95))
    
  # plot A: screening outcome (possible, impossible, present) across mutations
  p_mutation_screen <- ggplot(data_for_plotting) +
    geom_bar(aes(x = mutation_name, fill = mutation_category), width = 1)  +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.1, 
                                     hjust=0.95,
                                     size=10), 
          axis.text.y = element_text(size=10),
          axis.title.x=element_text(size=10, face = "bold"),
          axis.title.y=element_text(size=10, face = "bold"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=10),
          legend.position = "top") +
    labs(x = "Amino acid substitution",
         y = "Number of species") +
    scale_y_continuous(breaks = seq(0, 20000, by = 1000), expand = c(0.01, 0)) +
    scale_fill_manual(values = cols, name = "Mutation possibility:") +
    labs(fill = "") 
  
  ggsave(filename = file_name, p_mutation_screen, width = 8, height = 6)
}


plot_evolvability_by_class <-function(filtered_output, 
                                       bacterial_taxonomy, 
                                       n_classes_to_plot = 20,
                                       file_name) {
  # phylum abbreviations:
  phylum_abbreviations = c(Actinomycetota = "A",
                           Bacillota = "Ba",
                           Bacteroidota = "Bc",
                           Campylobacterota = "Ca",
                           Cyanobacteriota = "Cy",
                           Deinococcota = "D",
                           Mycoplasmatota = "M",
                           Planctomycetota = "Pl",
                           Pseudomonadota = "Ps",
                           Spirochaetota = "S",
                           Thermodesulfobacteriota = "T",
                           Misc = "Misc")
  #average of possibility per mutation across classes
  plot_data <- filtered_output |>
    rename(genus = genus) |>
    group_by(species, genus, accession_numbers, mutation_name, n_possible) |>
    summarise(n = n(), .groups = 'drop') |>
    filter(n == 1) |>            # remove multi-copy species
    left_join(bacterial_taxonomy |> select(genus, class)) |>
    filter(!is.na(class)) |>
    mutate(class = fct_lump_n(class, n = n_classes_to_plot)) |>
    group_by(class, mutation_name) |>
    summarise(n_pos = mean(n_possible), .groups = 'drop') |>
    left_join(bacterial_taxonomy |> select(class, phylum) |> distinct()) |>
    mutate(phylum = ifelse(class == "Other", "Misc", phylum)) |>
    mutate(phylum = phylum_abbreviations[phylum]) |>
    mutate(class = factor(class), phylum = factor(phylum)) |>
    mutate(phylum = fct_relevel(phylum, "Misc", after = Inf)) |>
    mutate(mutation_name = str_replace(mutation_name, "_", ""))
  
  p <- ggplot(plot_data) +
    geom_tile(aes(x = mutation_name, 
                  y = interaction(reorder(class, dplyr::desc(class)), 
                                  reorder(phylum, dplyr::desc(phylum)),
                                  sep = "!"), 
                  fill = n_pos)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
          #axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          #legend.title = element_text(size = 6),
          #legend.text = element_text(size = 6),
          legend.position = "bottom") +
    scale_y_discrete(guide = guide_axis_nested(delim = "!"), name = "Phylum and class") +
    scale_fill_gradient(low = hsv(0, 0, 0.2), high = hsv(0, 0, 1), limits = c(0, 3)) +
    labs(x = "Amino acid substitution",
         y = "Class",
         fill = "Mean evolvability")
  ggsave(filename = file_name, p, width = 5, height = 6)
}




#' produces a figure showing resistance and evolvability across classes and genera
#' @param filtered_output a data frame providing screened and filtered gene sequences
#' @param file_name the path that the plot should be save in
#' @param min_frac_resistant filter that gives the minimum fraction of species within
#' a genus that need to be resistant for this genus to be included in figure A. 
#' @param min_genus_size minimum number of species within a genus for this genus to be included in the figure 
#' @param n_genera_to_plot number of genera to be included in figure A. 
#' (Top n_genera_to_plot ranked by fraction of resistant species.)
#' @param n_muts_to_plot Number of mutations to be distinguished in Figure A. 
#' (The remaining ones will be lumped together in category "Other".)
#'
#' @return A plot as a ggplot object in pdf format
#' @export
#'
#' @examples plot_classes_genera(filtered_output, "./plot/myFilename.pdf")
plot_classes_genera <- function(filtered_output,
                                bacterial_taxonomy,
                                file_name,
                                n_classes_to_plot = 20,
                                min_frac_resistant = 0.1, 
                                min_genus_size = 12, 
                                n_genera_to_plot = 30,
                                n_muts_to_plot = 6) {
  
  # preparing the data:
 
  # merging species with multiple gene copies:
  merged_filtered_output <- filtered_output |>
    group_by(species, genus, accession_numbers, mutation_name) |>
    summarise(mutation_category = ifelse(any(mutation_category == "present"),
                                         "present",
                                         ifelse(any(mutation_category == "possible"), "possible", "impossible")), 
              .groups = "drop")
  
  # number of predicted mutations present in each species:
  n_muts_per_species <- merged_filtered_output |>
    group_by(species, accession_numbers, genus) |>
    summarise(muts_present = sum(mutation_category == "present"), 
              n_possible = sum(mutation_category == "possible"),
              # n_possible = ifelse(any(mutation_category == "present"),
              #                     NA, sum(mutation_category == "possible")),
              .groups = "drop")
  
  # table retaining only the first mutation for each species:
  first_mut_per_species <- merged_filtered_output |>
    filter(mutation_category == "present") |>
    group_by(species) |>
    summarise(first_mut = mutation_name[1])
  
  # table of species and which resistance mutations they have (including "none" or "multiple"):
  species_with_muts <- n_muts_per_species |>
    left_join(first_mut_per_species, by = join_by(species)) |>
    mutate(category = ifelse(muts_present == 0L, "None", ifelse(muts_present > 1, "Multiple", first_mut))) |>
    select(species, genus, category, n_possible) |>
    left_join(bacterial_taxonomy, by = join_by(genus == genus))|>
    distinct() |>
    mutate(class = factor(class)) |>
    mutate(genus = factor(genus)) |>
    mutate(category = str_replace(category, "_", "")) |>
    mutate(category = fct_lump_n(category, n = (n_muts_to_plot + 2))) |>
    mutate(category = fct_relevel(category, "None")) |>
    mutate(category = fct_relevel(category, "Multiple", after = Inf)) |>
    mutate(category = fct_recode(category, " " = "None"))
  
  cols <- c(" " = rgb(0,0,0,0),
            "43N" = brewer.pal(9, name = "Pastel1")[1],
            "43R" = brewer.pal(9, name = "Pastel1")[2],
            "43T" = brewer.pal(9, name = "Pastel1")[3],
            "86C" = brewer.pal(9, name = "Pastel1")[4],
            "91L" = brewer.pal(9, name = "Pastel1")[6],
            "88E" = brewer.pal(9, name = "Pastel1")[7],
            "88R" = brewer.pal(9, name = "Pastel1")[5],
            "Other" = "grey", 
            "Multiple" = "black")
  
  # Plot A: evolvability by class
  classes_for_plotting <- species_with_muts |>
    group_by(class) |>
    summarise(n = n(), .groups = "drop") |>
    filter(!is.na(class)) |>
    slice_max(n, n = n_classes_to_plot) |>
    pull(class)
  
  # plot_A1 <- ggplot(filter(species_with_muts, class %in% classes_for_plotting)) +
  #   geom_boxplot(aes(x = reorder(class, dplyr::desc(class)), 
  #                   y = n_possible), outliers = FALSE) +
  #   theme_bw() +
  #   scale_y_continuous(expand = c(0.01, 0)) +
  #   scale_fill_manual(values = cols, name = "") +
  #   labs(x = "Class",
  #        y = "Evolvability") +
  #   coord_flip()

  plot_A1 <- ggplot(filter(species_with_muts, class %in% classes_for_plotting)) +
  geom_violin(aes(x = reorder(class, dplyr::desc(class)), y = n_possible),
    width=1.5, size=0.3, alpha = 0.5
  ) + 
  theme_minimal() +
  scale_y_continuous(expand = c(0.01, 0)) +
    scale_fill_manual(values = cols, name = "") +
  labs(x = "Class", y = "Evolvability") +
  coord_flip()

  
  # Plot A2: predicted resistance mutations by class
  plot_A2 <- ggplot(filter(species_with_muts, class %in% classes_for_plotting)) +
    geom_bar(aes(x = reorder(class, dplyr::desc(class)), 
                 fill = category),
             position = "fill") +
    theme_bw() +
    scale_y_continuous(expand = c(0.01, 0), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = cols, name = "") +
    labs(x = "Class",
         y = "Resistant species") +
    coord_flip() +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) + 
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = margin(0.2, 0.4, 0.2, 0.2, "cm")
    )
      
  # Plot B: predicted resistance mutations by genus
  
  # genera to be plotted:
  genera_for_plotting <- species_with_muts |>
    group_by(genus) |>
    summarise(n = n(), n_res = sum(category != " "), .groups = "drop") |>
    filter(!(genus %in% c("Candidatus", "Wolbachia"))) |>
    mutate(fraction_res = n_res / n) |>
    filter(fraction_res >= min_frac_resistant, n >= min_genus_size) |>
    slice_max(fraction_res, n = n_genera_to_plot) |>
    pull(genus)
  
  plot_B <- ggplot(filter(species_with_muts, 
                          genus %in% genera_for_plotting) |>
                     mutate(genus = factor(paste0(" ", genus))) |>
                     mutate(class = as.character(class)) |>
                     mutate(class = ifelse(is.na(class), "?", substr(class, 1, 1))) |>
                     mutate(class = factor(class))) +
    geom_bar(aes(x = interaction(reorder(genus, dplyr::desc(genus)),
                                 reorder(class, dplyr::desc(class)), 
                                 sep = "!"), 
                 fill = category),
             position = "fill") +
    theme_bw() +
    scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Class and genus") +
    scale_y_continuous(expand = c(0.01, 0), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = cols, name = "") +
    labs(x = "Genus",
         y = "Resistant species") +
    coord_flip() +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    theme(plot.margin = margin(0.2, 0.3, 0.2, 0.7, "cm"))
  
  combined_plot <- ggarrange(plot_A1, plot_A2, plot_B,
                             ncol = 3,
                             widths = c(1, 0.55, 1),
                             labels = list("A", "", "  B"),
                             common.legend = TRUE,
                             legend = "bottom") + 
    theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
  
  ggsave(filename = file_name, combined_plot, width = 9, height = 7)
}


#' produces plots to show the distribution of gene sequence length, aligning score and distance from RRDR of ref seq
#'
#' @param final_output a data frame providing the results of mutation screening in all extracted gene sequences
#' @param filtered_output a data frame providing screened and filtered gene sequences
#' @param file_names the path that the plots should be save in
#'
#' @return A plot as a ggplot object in pdf format
#' @export
#'
#' @examples get_hist(final_output, target_sequences, "./plot/myFilename.pdf")
plot_target_sequences_stats <- function(final_output, 
                                        filtered_output, 
                                        min_seq_length,
                                        min_alig_score,
                                        max_core_dist,
                                        file_names) {
  cols <- c("raw" = "grey20", "filtered" = "grey80")
  
  dat_raw <- select(final_output, species, accession_numbers, gene_copy, target_length, alig_score, core_dist) |>
    mutate(filter = factor("raw")) |>
    distinct()
  dat_filtered <- select(filtered_output, species, accession_numbers, gene_copy, target_length, alig_score, core_dist) |>
    mutate(filter = factor("filtered")) |>
    distinct()
  dat <- rbind(dat_raw, dat_filtered)
  
  plot1 <- ggplot(dat) +
    geom_histogram(aes(target_length, fill = filter), position = "identity") +
    geom_vline(aes(xintercept = min_seq_length), col = "red") +
    theme_classic() +  # Apply a minimal theme
    theme(panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +
    labs(x = "Length", y = "Number of sequences", fill = "") +
    scale_fill_manual(values = cols)
  
  plot2 <- ggplot(dat) +
    geom_histogram(aes(alig_score, fill = filter), position = "identity", binwidth = 10) +
    geom_vline(aes(xintercept = min_alig_score), col = "red") +
    theme_classic() +  # Apply a minimal theme
    theme(panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +
    labs(x = "Alignment score", y = NULL, fill = "") +
    scale_fill_manual(values = cols)
  
  plot3 <- ggplot(dat) +
    geom_histogram(aes(core_dist, fill = filter), position = "identity", binwidth = 1) +
    geom_vline(aes(xintercept = max_core_dist), col = "red") +
    theme_classic() +  # Apply a minimal theme
    theme(panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +
    labs(x = "Levenshtein distance", y = NULL, fill = "") +
    scale_fill_manual(values = cols)
    
  hists <- ggarrange(plot1, plot2, plot3, 
                     nrow = 1, 
                     labels = LETTERS[1:3], 
                     common.legend = TRUE, vjust = 0)
  ggsave(filename = file_names[1], hists, width = 10, height = 4)
  
  # pairs plot showing covariation between the three variables:
  pairs_plot <- ggpairs(dat, 
                        columns = c("target_length", "alig_score", "core_dist"), 
                        aes(color = filter, alpha = 0.01),
                        columnLabels = c("length", "alignment score", "core distance")) +
    theme_bw()
  ggsave(filename = file_names[2], pairs_plot, width = 12, height = 12)
}


plot_multiseq_stats <- function(multiseq_stats, file_name) {
  cols <- c(
  brewer.pal(9, "Set1"),            
  "black"                         
)
  p <- ggplot(multiseq_stats) +
    geom_point(aes(x = alig_score, 
                   y = core_dist + rnorm(nrow(multiseq_stats), 0, 0.05),
                   colour = fct_lump_min(genus, 3),
                   shape = fct_relevel(case_when(resistance_status == 0 ~ "none", 
                                             resistance_status == 1 ~ "one copy", 
                                             resistance_status == 2 ~ "both copies"), "both copies", after = Inf))) +
    scale_x_continuous(limits = c(0, 6000), ) +
    scale_color_manual(values = cols) +
    labs(x = "Alignment score", y = "Levenshtein distance", colour = "Genus", shape = "Resistance") +
    theme_bw()
  ggsave(filename = file_name, p, width = 8, height = 5)
}


plot_cons <- function(cons, pos, dist_type, n_plots = 6, file_name) {

  mean_dsts <- cons$means |>
    mutate(pos = 1:nrow(cons$means))
  mean_dsts <- bind_rows(
    mutate(mean_dsts, pos_shift = pos - 0.4999),
    mutate(mean_dsts, pos_shift = pos + 0.4999)
  ) |>
    select(pos, pos_shift, everything()) |>
    arrange(pos_shift)
  
  if (dist_type == "grantham") {
    mean_dsts <- mean_dsts |>
      mutate(toplot_Ecoli = grantham_Ecoli) |>
      mutate(toplot_rnd = grantham_rnd)
  } else if (dist_type == "hamming") {
    mean_dsts <- mean_dsts |>
      mutate(toplot_Ecoli = hamming_Ecoli) |>
      mutate(toplot_rnd = hamming_rnd)
  } else {
    stop("Unknown dist_type argument. (Must be 'hamming' or 'grantham')")
  }

  pos_ranges <- split_indices(nrow(mean_dsts)/2, n_plots)
  p <- list()
  max_d <- max(c(mean_dsts$toplot_Ecoli, mean_dsts$toplot_rnd), na.rm = TRUE)
  
  for(i in 1:n_plots) {
    dat <- mean_dsts |>
      filter(pos_shift >= pos_ranges$start[i] & pos_shift <= pos_ranges$end[i])
    mut_pos <- pos[(pos <= pos_ranges$end[i]) & (pos >= pos_ranges$start[i]) ]
    p[[i]] <- ggplot() +
      geom_vline(xintercept = mut_pos, col = "red", linewidth = 1.2, alpha = 0.5) +
      geom_line(data = dat, mapping = aes(x = pos_shift, y = toplot_Ecoli)) +
      geom_line(data = dat, mapping = aes(x = pos_shift, y = toplot_rnd), col = "blue") +
      xlab("Position") +
      ylab("AA diversity") +
      ylim(c(0, max_d)) +
      theme_bw()
  }
  
  p_combined <- p[[1]]
  for(i in 2:n_plots) {
    p_combined <- p_combined + p[[i]]
  }
  p_combined <- p_combined + plot_layout(ncol = 1)
  ggsave(filename = file_name, p_combined, width = 9, height = 12)
}



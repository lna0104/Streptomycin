# This file contains functions for producing plots and figures.


#' Produces figure giving the frequency of reported position, mutation and species based on Google sheet
#' @param muts a data frame table providing info on all reported mutations in the literature
#' @param file_name the path that the plot should be save in
#' @param n_frequency number of reported mutations, positions and species
#' @return Three plots as ggplot objects in one single pdf file
#' @export
#'
#' @examples plot_reported_mutations_nt(muts, "./plot/myFilename.pdf")
plot_reported_mutations_nt <- function(muts, file_name, n_frequency) {

  #1. the frequency of reported mutations across sites (currently not included in figure)
  frequency_position <-  muts %>% 
    select(Species, Nt_pos_Ecoli, Origin) %>%
    filter(!is.na(Nt_pos_Ecoli)) %>%
#    filter(AA_pos_Ecoli %in% mutation_list_reports$AA_pos_Ecoli) %>% 
    group_by(Species, Nt_pos_Ecoli) %>%
    summarise(Origin = ifelse(all(Origin == "Isolate"), "Isolate",
                              ifelse(all(Origin == "Lab-generated"), "Lab-generated", "Both")),
              .groups = "drop")

  plot1 <- frequency_position %>%
    mutate(Nt_pos_Ecoli = factor(Nt_pos_Ecoli, levels = sort(unique(Nt_pos_Ecoli)))) %>%
    ggplot() +
    geom_bar(aes(x = Nt_pos_Ecoli, fill = Origin), width = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.1, 
                                     size = 10),
          axis.text.y = element_text(size = 10), 
          axis.title=element_text(size=10, face="bold"),
          plot.title=element_text(hjust=0),
          legend.position = "top")+
    labs(x = "Nucleotide position",
         y = "Number of species", 
         fill = "Origin:") +
    scale_y_continuous(breaks = seq(0, 50, by = 5), expand = c(0.01, 0)) +
    scale_fill_brewer(palette = "Pastel1") 

  #2. the frequency of each reported amino acid substitution across species:

  frequency_mutation <- muts |>
    select(Species, Nt_pos_Ecoli, Nt_mutation, Origin) |>
    filter(!is.na(Nt_pos_Ecoli), !is.na(Nt_mutation)) |>
    mutate(mutation_name = paste0(Nt_pos_Ecoli, Nt_mutation)) |> 
    select(Species, mutation_name, Origin) %>%
    group_by(Species, mutation_name) %>%
    summarise(Origin = ifelse(all(Origin == "Isolate"), "Isolate",
                              ifelse(all(Origin == "Lab-generated"), "Lab-generated", "Both")),
              .groups = "drop") |>
    # arrange based on numeric position in mutation name
    mutate(pos_num = as.numeric(str_extract(mutation_name, "^[0-9]+"))) |>
    arrange(pos_num, mutation_name) |>
    mutate(mutation_name = factor(mutation_name, 
                                  levels = unique(mutation_name))) |>
    group_by(mutation_name) |>
    filter(n() >= n_frequency)
  
  plot2 <- ggplot(frequency_mutation) +
    geom_bar(aes(x = mutation_name, fill = Origin))  +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.1, 
                                     hjust=0.95, 
                                     size=8), 
          axis.text.y = element_text(size=8), 
          axis.title=element_text(size=10, face="bold"),
          plot.title=element_text(hjust=0),
          legend.title = element_text(size=10),
          legend.text = element_text(size=8),
          legend.position = "top") +
    labs(x = "Nucleotide substitution",
         y = "Number of species", 
         fill = "Origin:") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_fill_manual(values = c(
      "Lab-generated" = "#3B9AB2",  
      "Isolate" = "#EBCC2A",       
      "Both" = "#EF5703"            
    ), breaks=c('Lab-generated', 'Isolate', 'Both')) 

  #3. how many mutations for each species have been reported
  frequency_per_species <-  muts |>
    select(Species, Nt_pos_Ecoli, Nt_mutation, Origin) |>
    filter(!is.na(Nt_pos_Ecoli), !is.na(Nt_mutation)) |>
    mutate(mutation_name = paste0(Nt_pos_Ecoli, Nt_mutation)) |>
    group_by(mutation_name, Species) |>
    summarise(Origin = ifelse(all(Origin == "Isolate"), "Isolate",
                              ifelse(all(Origin == "Lab-generated"), "Lab-generated", "Both")),
              .groups = "drop") |>
    group_by(Species) |> 
    filter(n() >= n_frequency) 
  
  plot3 <- frequency_per_species %>%
    ggplot() +
    geom_bar(aes(x = fct_infreq(Species), fill = Origin), width = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     hjust=0.9,
                                     vjust=0.2, 
                                     size=8, face = "italic"), 
          axis.title.y=element_text(size=10, face="bold", hjust=1.2),
          axis.title.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10), 
          plot.title=element_text(hjust=0), 
          legend.position = "none") +
    labs(x = "Species",
         y = "Number of reported mutations", 
         fill = "Origin:") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_fill_manual(values = c(
      "Lab-generated" = "#3B9AB2",  
      "Isolate" = "#EBCC2A",       
      "Both" = "#EF5703"            
    ), breaks=c('Lab-generated', 'Isolate', 'Both')) 
  
  # ggsave(filename = "./plots/rrs_reported_mutations_per_species.pdf", plot3, width = 10, height = 6)

  #combine plots:
  reported_mutations <- ggarrange(plot2, plot3, 
                                  nrow=2, common.legend = FALSE,
                                  labels = c("A", "B"), hjust = 0, vjust = 1.2)
  
  ggsave(filename = file_name, reported_mutations, width = 6, height = 6)
}


plot_reported_article_before_process <- function(muts, n_frequency = 1, file_name){
  plotting_data <- muts |>
    mutate(Mutation_name = paste0(Nt_original, Nt_pos, Nt_mutation)) |>
    select(Species, Mutation_name, Ref_code) |>
    distinct()
  
  n_articles_per_mutation <- plotting_data |>
    group_by(Mutation_name) |>
    summarise(n_articles = n_distinct(Ref_code)) |>
    filter(n_articles>n_frequency) |>
    arrange(n_articles)


  plot1 <- n_articles_per_mutation %>%
    ggplot(aes(x = reorder(Mutation_name, -n_articles), y = n_articles)) +
    geom_bar(stat = "identity", fill = "grey75") +
    labs(
      x = "Mutations",
      y = "Number of articles"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.1, 
                                     hjust=0.95, 
                                     size=8), 
          axis.text.y = element_text(size=8), 
          axis.title=element_text(size=10, face="bold"),
          plot.title=element_text(hjust=0),
          legend.title = element_text(size=10),
          legend.text = element_text(size=8)) 
  
  n_articles_per_species <- plotting_data |>
    group_by(Species) |>
    summarise(n_articles = n_distinct(Ref_code)) |>
    # filter(n_articles>n_frequency) |>
    arrange(n_articles)
  
  plot2 <- n_articles_per_species %>%
    ggplot(aes(x = reorder(Species, -n_articles), y = n_articles)) +
    geom_bar(stat = "identity", fill = "grey75") +
    labs(
      x = "Species",
      y = "Number of articles"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.1, 
                                     hjust=0.95, 
                                     size=8), 
          axis.text.y = element_text(size=8), 
          axis.title=element_text(size=10, face="bold"),
          plot.title=element_text(hjust=0),
          legend.title = element_text(size=10),
          legend.text = element_text(size=8)) 
  
  #combine plots:
  combine_plots <- ggarrange(plot1, plot2, 
                                  nrow=2, heights = c(0.8, 1.2),
                                  common.legend = FALSE,
                                  labels = c("A", "B"))
  
  ggsave(filename = file_name, combine_plots, width = 12, height = 8)

  
}

#' Produces figure giving the frequency of mutations in present, possible and impossible categories
#' @param filtered_output  a data frame providing screened and filtered gene sequences
#' @param file_name the path that the plot should be save in
#' @return A plot as a ggplot object in pdf format
#' @export
#'
#' @examples mutation_screen(filtered_output, "myFilename.pdf")
plot_mutation_screen_nt <- function(filtered_output, file_name) {
  
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
    # mutate(mutation_category = ifelse(mutation_category == "possible", 
    #                                   n_possible, 
    #                                   mutation_category)) |>
    mutate(mutation_category = factor(mutation_category, 
                                      levels = c("possible", "present"))) |>
    mutate(mutation_name = str_replace(mutation_name, "_", "")) |>
    mutate(mutation_name = factor(mutation_name, levels = sort(unique(mutation_name))))
  
  # cols = c(hsv(0, 0, c(0.2, 0.73, 0.8)), hsv(0, 1, 0.95))
  # cols <- c("#899DA4","#f9e5ae", "#FAEFD1", "#C93312")
  cols <- c("#f9e5ae", "#C93312")
  
  # plot A: screening outcome (possible, impossible, present) across mutations
  p_mutation_screen <- ggplot(data_for_plotting) +
    geom_bar(aes(x = mutation_name, fill = mutation_category), width = 1, 
             color = "white", linewidth=0.2)  +
    theme_bw() +
    theme(axis.text.x = element_text(angle=0, 
                                     vjust=0.5, 
                                     hjust=0.95,
                                     size=10), 
          axis.text.y = element_text(size=10),
          axis.title.x=element_text(size=10, face = "bold"),
          axis.title.y=element_text(size=10, face = "bold"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=10),
          legend.position = "bottom") +
    labs(x = "Amino acid substitution",
         y = "Number of species") +
    scale_y_continuous(breaks = seq(0, 20000, by = 1000), expand = c(0.01, 0)) +
    scale_fill_manual(values = cols, name = "Mutation possibility:") +
    labs(fill = "") + 
    coord_flip()
  
  ggsave(filename = file_name, p_mutation_screen, width = 8, height = 4)
}

#' produces a figure showing resistance and evolvability across classes and genera
#' @param filtered_output a data frame providing screened and filtered rpoB sequences
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
#' @examples plot_classes_genera_nt(filtered_output, "./plot/myFilename.pdf")
plot_classes_genera_nt <- function(filtered_output,
                                bacterial_taxonomy,
                                file_name,
                                n_classes_to_plot = 20,
                                #min_frac_resistant = 0.1, 
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
  
  # restrict to species that actually have a mutation 
  species_with_any_mut <- species_with_muts |>
    dplyr::filter(category != " ")
  
  cols <- c(" " = rgb(0,0,0,0),
            "523C" = wes_palette("FantasticFox1")[3],
            "885A" = wes_palette("FantasticFox1")[1],
            "912G" = wes_palette("FantasticFox1")[2],
            "913G" = wes_palette("AsteroidCity2")[3],
            "Other" = "grey", 
            "Multiple" = wes_palette("FantasticFox1")[5])
  
  # classes_for_plotting <- species_with_muts |>
  #   group_by(class) |>
  #   summarise(n = n(), .groups = "drop") |>
  #   filter(!is.na(class)) |>
  #   slice_max(n, n = n_classes_to_plot) |>
  #   pull(class)
  
  classes_for_plotting <- species_with_any_mut |>
    group_by(class) |>
    summarise(n = n(), .groups = "drop") |>
    filter(!is.na(class)) |>
    slice_max(n, n = n_classes_to_plot) |>
    pull(class)
  
  # Plot: predicted resistance mutations by class
  bar_plot <- ggplot(filter(species_with_muts, class %in% classes_for_plotting)) +
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
      legend.position = "bottom",
      plot.margin = margin(0.2, 0.4, 0.2, 0.2, "cm")
    )

  ggsave(filename = file_name, bar_plot, width = 16, height = 6)
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
plot_target_sequences_stats_nt <- function(final_output, 
                                        filtered_output, 
                                        min_seq_length,
                                        min_alig_score,
                                        max_core_dist,
                                        file_names) {
  big_txt <- theme(
    axis.title = element_text(size = 14, face="bold"),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 13)
  )
  
  cols <- c("raw" = "grey20", "filtered" = "grey80")
  
  dat_raw <- select(final_output, species, accession_numbers, gene_copy, target_length, alig_score, core_dist) |>
    mutate(filter = factor("raw")) |>
    distinct()
  dat_filtered <- select(filtered_output, species, accession_numbers, gene_copy, target_length, alig_score, core_dist) |>
    mutate(filter = factor("filtered")) |>
    distinct()
  dat <- rbind(dat_raw, dat_filtered)
  
  plot1 <- ggplot(dat) +
    geom_histogram(aes(target_length, fill = filter), position = "identity",
                   binwidth = 100) +
    geom_vline(aes(xintercept = min_seq_length), col = "red") +
    theme_classic() +  # Apply a minimal theme
    theme(panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +
    big_txt +
    labs(x = "Length", y = "Number of sequences", fill = "") +
    scale_fill_manual(values = cols)
  
  plot2 <- ggplot(dat) +
    geom_histogram(aes(alig_score, fill = filter), position = "identity", binwidth = 10) +
    geom_vline(aes(xintercept = min_alig_score), col = "red") +
    theme_classic() +  # Apply a minimal theme
    theme(panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +
    big_txt +
    labs(x = "Alignment score", y = NULL, fill = "") +
    scale_fill_manual(values = cols)
  
  plot3 <- ggplot(dat) +
    geom_histogram(aes(core_dist, fill = filter), position = "identity", binwidth = 1) +
    geom_vline(aes(xintercept = max_core_dist), col = "red") +
    theme_classic() +  # Apply a minimal theme
    theme(panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +
    big_txt +
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


plot_multiseq_stats_nt <- function(multiseq_stats, file_name) {
  p <- ggplot(multiseq_stats) +
    geom_point(aes(x = mean_alig_score, 
                   y = mean_core_dist), color = "steelblue") +
    scale_x_continuous(limits = c(0, 5000), ) +
    labs(x = "Mean alignment score", y = "Mean Levenshtein distance") +
    theme_bw()
  ggsave(filename = file_name, p, width = 8, height = 5)
}



plot_multiseq_dist <- function(filtered_output, file_name){
  
  plotting_data <- filtered_output |>
    distinct(accession_numbers, target_name, gene_copy) |>
    count(accession_numbers, name = "n")
  
  dist_plot <- ggplot(plotting_data, 
                      aes(x=n)) + 
    geom_histogram(
      binwidth = 1,
      boundary = 0,
      colour = "grey50", 
      fill = "grey80",
    ) +
    scale_x_continuous(
      limits = c(1, NA),  
      expand = c(0, 1), 
      breaks = function(lims) seq(floor(lims[1]), ceiling(lims[2]), by = 1) + 0.5,
      labels = function(x) as.character(floor(x))
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14, face="bold"),
      axis.text  = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 13)    
      ) +
    labs(
      x = "Gene copies",
      y = "Number of species"
    )
    
  ggsave(filename = file_name, width = 10, height = 8)
}




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
                                     size=8), 
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
#' @param filtered_output a data frame providing screened and filtered gene sequences
#' @param file_name the path that the plot should be save in
#' @param genus_variants a data frame identifying genera in the filtered_output that correspond to multiple entries in the GTDB taxonomy data (e.g., "Actinomadura" represented as "Actinomadura_C" and "Actinomadura_D").
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
                                gtdb_taxonomy,
                                genus_variants,
                                file_name,
                                n_classes_to_plot = 20,
                                min_frac_resistant = 0.1, 
                                min_genus_size = 10, 
                                n_genera_to_plot = 20,
                                n_muts_to_plot = 6) {
  
  # preparing the data:
  
  # merging species with multiple gene copies:
  merged_filtered_output <- filtered_output |>
    group_by(species, genus, accession_numbers, mutation_name) |>
    summarise(
        mutation_category = if (any(mutation_category == "present")) "present" else "possible",
        .groups = "drop"
      )
  
  # number of predicted mutations present in each species:
  n_muts_per_species <- merged_filtered_output |>
    group_by(species, genus, accession_numbers) |>
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
  processed_data <-  n_muts_per_species |>
    left_join(first_mut_per_species, by = join_by(species)) |>
    mutate(category = ifelse(muts_present == 0L, "None", ifelse(muts_present > 1, "Multiple", first_mut))) |>
    select(species, genus, category, n_possible) |>
    left_join(gtdb_taxonomy, by = join_by(genus == genus), relationship = "many-to-many" )|>
    left_join(genus_variants |> select(genus_origin, class_var = class, phylum_var = phylum),
              by = join_by(genus == genus_origin), relationship = "many-to-many") |> 
    mutate(
      class = if_else(is.na(class), class_var, class),
      phylum = if_else(is.na(phylum), phylum_var, phylum)
    ) |>
    select(-class_var, -phylum_var) |>
    distinct() 
  
  species_with_muts <- processed_data |>
    mutate(class = factor(class)) |>
    mutate(genus = factor(genus)) |>
    mutate(category = str_replace(category, "_", "")) |>
    mutate(category = fct_lump_n(category, n = (n_muts_to_plot + 2))) |>
    mutate(category = fct_relevel(category, "None")) |>
    mutate(category = fct_relevel(category, "Multiple", after = Inf)) |>
    mutate(category = fct_recode(category, " " = "None"))
  
  # cols <- c(" " = rgb(0,0,0,0),
  #           "43N" = brewer.pal(9, name = "Pastel1")[1],
  #           "43R" = brewer.pal(9, name = "Pastel1")[2],
  #           "43T" = brewer.pal(9, name = "Pastel1")[3],
  #           "86C" = brewer.pal(9, name = "Pastel1")[4],
  #           "91L" = brewer.pal(9, name = "Pastel1")[6],
  #           "88E" = brewer.pal(9, name = "Pastel1")[7],
  #           "88R" = brewer.pal(9, name = "Pastel1")[5],
  #           "Other" = "grey", 
  #           "Multiple" = "black")
  
  cols <- c(" " = rgb(0,0,0,0),
            "523C" =  wes_palette("Zissou1")[5],
            "885A" =  wes_palette("Zissou1")[3],
            "912G" =  wes_palette("Zissou1")[1],
            "913G" =  wes_palette("Zissou1")[2],
            "Other" = "grey",
            "Multiple" = "black")
  
  # Plot A: Number species per class
  n_species_per_class <- species_with_muts |>
    mutate(class = as.character(class), 
           class = ifelse(is.na(class) | class == "", "Unidentified", class)) |>
    group_by(class) |>
    summarise(n = n(), .groups = "drop")
  
  classes_for_plotting <- n_species_per_class |>
    filter(class != "Unidentified") |>
    slice_max(n, n = n_classes_to_plot) |>
    pull(class)
  
  plot_A1 <- ggplot(filter(species_with_muts, class %in% classes_for_plotting)) +
    # geom_jitter(
    # aes(x = reorder(class, dplyr::desc(class)), y = n_possible),
    # width = 0.2, size = 2.5, alpha = 0.1
    # ) +
    geom_violin(aes(x = reorder(class, dplyr::desc(class)), y = n_possible),
                width=1.2, size=0.3
    ) + 
    # scale_fill_manual(values = cols_class, guide = "none") + 
    scale_y_continuous(expand = c(0.01, 0)) +
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
    labs(x = "Class",
         y = "Resistant species") +
    coord_flip() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank()) + 
    scale_fill_manual(values = cols, name = " ") +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) 
  
  plot_A2_no_legend <- plot_A2 + theme(legend.position = "none")
  legend_a2 <- get_legend(plot_A2)
  
  # Plot B: predicted resistance mutations by genus
  
  # genera to be plotted:
  genera_for_plotting <- species_with_muts |>
    group_by(genus) |>
    summarise(n = n(), n_res = sum(category != " "), .groups = "drop") |>
    # filter(!(genus %in% c("Candidatus", "Wolbachia"))) |>
    mutate(fraction_res = n_res / n) |>
    filter(fraction_res >= min_frac_resistant, n >= min_genus_size) |>
    slice_max(fraction_res, n = n_genera_to_plot) |>
    pull(genus)
  
  # Prepare taxonomy levels
  taxonomy <- species_with_muts |> 
    filter(genus %in% genera_for_plotting) |> 
    select(-c(category, n_possible)) 
  
  # Generate edges from class → order → family → genus 
  taxonomy_edges <- bind_rows(
    taxonomy |> transmute(from = class, to = order),
    taxonomy |> transmute(from = order, to = family),
    taxonomy |> transmute(from = family, to = genus)
  )|> 
    distinct() |>
    filter(!is.na(from), !is.na(to))
  
  # Create graph object
  graph <- tbl_graph(edges = taxonomy_edges, directed = TRUE)
  
  # # Function to propagate class name down the tree
  # propagate_class <- function(graph_tbl, class_names) {
  #   V(graph_tbl)$class_parent <- NA
  
  #   for (class_node in which(V(graph_tbl)$name %in% class_names)) {
  #     descendants <- igraph::subcomponent(graph_tbl, class_node, mode = "out")
  #     V(graph_tbl)$class_parent[descendants] <- V(graph_tbl)$name[class_node]
  #   }
  
  #   graph_tbl
  # }
  
  # Apply propagation
  # class_names <- names(cols_class)
  # graph_colored <- propagate_class(graph, class_names)
  
  # Annotate class and genus in the plot
  graph_for_plotting <- graph |> 
    activate(nodes) |> 
    mutate(
      rank = case_when(
        name %in% unique(taxonomy$class) ~ "class",
        name %in% unique(taxonomy$genus) ~ "genus",
        TRUE ~ "Other"
      )) |>
    # color = cols_class[class_parent]) |>
    arrange(desc(rank == "class"), desc(name))
  
  # Plot the multi-layer tree
  plot_B1 <- ggraph(graph_for_plotting, layout = "sugiyama") +
    geom_edge_link() +
    geom_node_point(aes(filter = rank == "Other"), alpha = 0) +  
    geom_node_label(aes(label = ifelse(rank %in% c("class", "genus"), name, "")), repel = TRUE) +
    theme_void() +
    # scale_colour_manual(values = cols_class) +
    scale_y_reverse() +
    coord_flip() +
    theme(
      plot.margin = margin(0, 0.3, 0, 0.2, "cm")
    )
  
  
  #Extract genus order
  layout <- create_layout(graph_for_plotting, layout = "sugiyama")

  genus_order <- as_tibble(layout, active = "nodes") |> 
    filter(rank == "genus") |>
    arrange(x) |> 
    pull(name) 
  
  plot_B2 <- ggplot(filter(species_with_muts, 
                           genus %in% genus_order) |>
                      mutate(genus = factor(genus, levels=genus_order)) |>
                      mutate(class = as.character(class)) |>
                      #  mutate(class = ifelse(is.na(class), "?", substr(class, 1, 1))) |>
                      mutate(class = factor(class))) +
    geom_bar(aes(x = genus, fill = category),
             position = "fill") +
    theme_bw() +
    # scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Class and genus") +
    scale_y_continuous(expand = c(0.01, 0), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = cols, name = "") +
    labs(y = "Resistant species") +
    coord_flip() +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    theme(
      plot.margin = margin(0.5, 0.5, 0.2, 0, "cm"), 
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
  
  # Combine A1 and A2 
  plot_A <- plot_grid(plot_A1, plot_A2_no_legend, 
                      nrow = 1, rel_widths = c(1.2, 1))
  
  #Combine B1 and  B2 vertical
  plot_B <- plot_grid(plot_B1, plot_B2,
                      ncol = 2,
                      rel_widths = c(1.2, 1),
                      labels = c("", ""),
                      align = "h") +
    theme(
      plot.margin = margin(0, 0.2, 0, 0.5, "cm")
    )         
  
  
  main_plot <- plot_grid(plot_A, plot_B,
                         ncol = 2,
                         rel_widths = c(1, 1),
                         labels = c("", "C"))
  
  # Add the legend 
  add_legend_plot <- plot_grid(main_plot, legend_a2, ncol = 1, rel_heights = c(1, 0.07)) 
  
  # combined_plot <- ggarrange(plot_B1, plot_B2, plot_C,
  #                            ncol = 3,
  #                            widths = c(1, 0.6, 1),
  #                            labels = list("A", "", "  B"),
  #                            common.legend = TRUE,
  #                            legend = "bottom") + 
  #   theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
  
  ggsave(filename = file_name, add_legend_plot, width = 14, height = 8)
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


# Plot number of species across different gene copy counts
# Species are colored by whether mutations are present 
plot_mutation_screen_gene_copies <- function(filtered_output, file_name){

  gene_copies_species <- filtered_output |>
    select(accession_numbers, species, target_name, mutation_name, mutation_category) |>
    distinct() |>
    group_by(accession_numbers, species, target_name) |>
    summarise(
      mutation_category_per_mut = if (any(mutation_category == "present")) {
      "present"
    } else {
      "not present"
    },
    .groups = "drop"
    ) |>
    group_by(accession_numbers, species) |>
    summarise(
      gene_copy = n(),
      mutation_category_per_species = if (any(mutation_category_per_mut == "present")) {
        "present"
      } else {
        "not present"
      },
      .groups = "drop"
    )


  # plot histogram number of species with and with out mutations per gene copies 
  hist_plot <- ggplot(gene_copies_species, aes(x = gene_copy, fill = mutation_category_per_species)) +
    geom_histogram(binwidth = 1) +
    scale_x_continuous(
      limits = c(0, 37),
      breaks = 1:37, 
      expand = c(0, 0.1)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(
      values = c(
        "present" = "#C93312",   
        "not present"   = "#899DA4"    
      )
    ) +
    labs(
        x = "Gene copies",
        y = "Number of genomes",
        fill = "Mutation category") +
    theme_classic()+
    theme(
      axis.title=element_text(size=10, face="bold"),
      axis.text = element_text(size = 8)
    )
  
  ggsave(filename = file_name, hist_plot, width = 12, height = 8)
}

# Plot gene copy distribution per mutation
# Species are colored by whether mutations are present 

plot_gene_copies_per_mutation <- function(filtered_output, file_name){
  
  n_copy_per_species <- filtered_output |>
    select(accession_numbers, species, target_name) |>
    distinct() |>
    group_by(accession_numbers, species) |>
    summarise(
      gene_copy = n(),
      .groups = "drop"
    ) 
  
  gene_copies_per_mut <- filtered_output |>
    group_by(accession_numbers, species, mutation_name) |>
    summarise(
      mutation_category_per_species = if (any(mutation_category == "present")) {
        "present"
      } else {
        "not present"
      },
      .groups = "drop"
    ) |> 
    left_join(n_copy_per_species, by = c("accession_numbers", "species")) 
  
  hist_plot <- ggplot(gene_copies_per_mut, aes(x = gene_copy, fill = mutation_category_per_species)) +
    geom_histogram(binwidth = 1) +
    scale_x_continuous(
      limits = c(0, 37),
      breaks = 1:37, 
      expand = c(0, 0.1)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(
      values = c(
        "present" = "#C93312",   
        "not present" = "#899DA4"    
      )
    ) +
    labs(
      x = "Gene copies",
      y = "Number of genomes",
      fill = "Mutation category"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size=10, face="bold"),
      axis.text  = element_text(size=6),
      strip.text = element_text(size=9, face="bold")
    ) +
    facet_wrap(~ mutation_name, scales = "free_y")
  
  ggsave(filename = file_name, hist_plot, width = 12, height = 10)
}

plot_variance_gene_copies <- function(filtered_output, file_name){
  n_copy_per_species <- filtered_output |>
    select(accession_numbers, species, target_name) |>
    distinct() |>
    group_by(accession_numbers, species) |>
    summarise(
      gene_copy = n(),
      .groups = "drop"
    ) |>
    filter(gene_copy > 1) 
  
  n_variants_per_species <- filtered_output |>
    filter(accession_numbers %in% n_copy_per_species$accession_numbers) |> 
    group_by(accession_numbers, species, Nt_pos_Ecoli) |> 
    summarise(
      n_variants = n_distinct(Nt_target),
      .groups = "drop"
    ) 
  
  bar_plot <- ggplot(n_variants_per_species, aes(x = factor(n_variants))) +
    geom_bar(fill = "lightgray", color = "black") +
    facet_wrap(~ Nt_pos_Ecoli, ncol = 2) + 
    labs(
      x = "Number of unique nucleotides",
      y = "Number of genomes"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 8)
    )
  
  ggsave(filename = file_name, bar_plot, width = 6, height = 6)

}

#' Produce bar plots showing species counts or percentages per class and phylum
#'
#' This function generates bar plots of species distribution across taxonomic 
#' classes and phyla, grouped by a categorical variable (e.g., "rrs", "rpsL", "both").
#' It supports nested x-axis labels for class/phylum and can plot either absolute 
#' counts or percentages of species.
plot_class_counts_and_percents <- function(processed_data, file_name){
  
  # Make a combined dataset: raw counts + proportions
  plot_data_counts <- processed_data |>
    mutate(plot_type = "Count")
  
  plot_data_props <- processed_data |>
    group_by(class, phylum) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    mutate(plot_type = "Proportion")
  
  # Stack them together
  plot_data <- bind_rows(
    plot_data_counts |> mutate(value = n),
    plot_data_props |> mutate(value = prop)
  )
  
  ggplot(plot_data, aes(
    x = interaction(class, phylum, sep = "!"),
    y = value,
    fill = category
  )) +
    geom_col() +
    scale_x_discrete(guide = guide_axis_nested(key = "!"), name = "Phylum and class") +
    scale_fill_manual(values = c(
      "rrs" = '#21908dff',
      "rpsL" = '#fde725ff',
      "both" = '#440154ff'
    )) +
    facet_grid(rows = vars(plot_type), scales = "free_y") +
    labs(
      x = "Class",
      fill = "Category"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      ggh4x.axis.nestline.x = element_line(linetype = 2),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = -0.01, hjust = 1, size = 8),
      axis.title.x = element_text(size = 10)
    )
  ggsave(filename = file_name, width = 14, height = 8)
}




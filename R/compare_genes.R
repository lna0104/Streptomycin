library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(ggh4x)
library(legendry)



# Load data
rrs_gene_list <- read_csv("./results/rrs_predicted_resistance.csv")
rpsL_gene_list <- read_csv("./results/predicted_resistance.csv")

# Combine all species 
all_species <- bind_rows(rrs_gene_list, rpsL_gene_list) %>%
  distinct(species, genus)

# Venn diagram showing the overlapping number of resistant species
venn.diagram(
  x = list(rrs_gene_list |> pull(species), 
           rpsL_gene_list |> pull(species)),
  category.names = c("rrs", "rpsL"),
  filename = './plots/venn_diagram.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 550 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c('#21908dff', '#fde725ff'),
  fill = c(alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.8,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c('#21908dff', '#fde725ff')
)

# Classify category
combined_species <- all_species %>%
  mutate(
    category = case_when(
      species %in% rrs_gene_list$species & species %in% rpsL_gene_list$species ~ "both",
      species %in% rrs_gene_list$species ~ "rrs",
      species %in% rpsL_gene_list$species ~ "rpsL"
    )
  )

# Generate gtdb_taxonomy file from GTDB metadata
gtdb_taxonomy <- read_csv("./data/gtdb_taxonomy.csv", show_col_types = FALSE) #bacterial taxonomic information from gtdb
genus_variants <- read_csv("./output/variants_gtdb_taxonomy.csv", show_col_types = FALSE)

processed_data <- combined_species |>
  select(species, genus, category) |>
  left_join(gtdb_taxonomy, by = join_by(genus == genus))|>
  left_join(genus_variants |> select(genus_origin, class_var = class, phylum_var = phylum),
            by = join_by(genus == genus_origin), relationship = "many-to-many") |> 
  mutate(
    class = if_else(is.na(class), class_var, class),
    phylum = if_else(is.na(phylum), phylum_var, phylum)
  ) 

p<-ggplot(processed_data, aes(
      x = interaction(class,
                      phylum,
                      sep = "!"),  
      fill = category
    )) +
      geom_bar() +
      guides(x = legendry::guide_axis_nested(key = "!")) +
      scale_fill_manual(values = c(
        "rrs" = '#21908dff',
        "rpsL" = '#fde725ff',
        "both" = '#440154ff'
      )) +
      labs(
        x = "Phylum and Class",
        y = "Number of species",
        fill = "Category"
      ) +
      theme_bw() +
      theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = -0.01, hjust = 1, size = 8),
        axis.title = element_text(size = 10)
    )

ggsave(
  filename = "./plots/phylum_class_rrs_rpsL.pdf",
  plot = p,
  width = 10,   
  height = 6
)




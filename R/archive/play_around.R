#sum of average possibility per class
num_poss <- filtered_output %>%
  rename(genus = genus_name) %>%
  full_join(bacterial_taxonomy) %>%
  replace_na(list(family = 'other',
                  order = 'other',
                  class = "other",
                  phylum = "other")) %>%
  filter(!is.na(AA_pos_Ecoli)) %>%
  group_by(class, mutation_name) %>%
  summarise(n_pos = mean(n_possible), .groups = 'drop') %>%
  group_by(class) %>%
  summarise(t_pos = sum(n_pos)) %>%
  arrange(desc(t_pos))


t_mutation_possibility <- ggplot(num_poss) +
  geom_line(aes(x = reorder(class, -t_pos), y = t_pos, group = 1), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 0.95, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  labs(x = "Class",
       y = "Sum of average possibility")+
  scale_y_continuous(breaks = seq(30, 60, by = 2), expand = c(0.01, 0))
t_mutation_possibility

#average of possibility per mutation across classes
num_poss <- filtered_output %>%
  rename(genus = genus_name) %>%
  full_join(bacterial_taxonomy) %>%
  replace_na(list(family = 'other',
                  order = 'other',
                  class = "other",
                  phylum = "other")) %>%
  filter(!is.na(AA_pos_Ecoli)) %>%
  group_by(class, mutation_name) %>%
  summarise(n_pos = mean(n_possible), .groups = 'drop') 

p_mutation_possibility <- ggplot(num_poss, aes(x = mutation_name, y = class)) +
  geom_tile(aes(fill = n_pos)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "top") +
  labs(x = "Amino acid substitution",
       y = "Class",
       fill = "Average of possibility")+
  scale_fill_gradient(low = "white", high = "dodgerblue4")

p_mutation_possibility
#####################################
#mutation presence per class
mut_present <- filtered_output %>% 
  rename(genus = genus_name) %>% 
  full_join(bacterial_taxonomy) %>% 
  replace_na(list(family = 'other', 
                  order = 'other',
                  class = "other",
                  phylum = "other")) %>% 
  filter(!is.na(AA_pos_Ecoli)) %>% 
  filter(mutation_category == "present") %>% 
  group_by(class, mutation_name) %>% 
  summarise(mut_pres = n(), .groups = 'drop') %>% 
  filter(mut_pres > 6)

p_mutation_presence <- ggplot(mut_present) +
  geom_bar(aes(x = class, y = mut_pres, fill = mutation_name),
           stat = "identity", width = 1, position = "fill")  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.1,
                                   hjust = 0.95,
                                   size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold")) +
  labs(x = "Class",
       y = "Fraction of mutation presence",
       fill = "Mutation:") +
  scale_y_continuous(labels = scales::number_format(scale = 1)) +
  scale_fill_brewer(palette = "Set3")

p_mutation_presence
#######Import data 
install.packages("tidyverse")
library(ggplot2)
Abund_strong<-read.table("data/bwa_counts_total_filtered_wMetadata_15mil_s95_strongFilter_relAbund.tsv", sep = '\t', header = TRUE)

nrow(Abund_strong) #1598 mags
library(tidyverse)
Abund_strong_2 <- Abund_strong %>% 
  separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(across(everything(), ~ sub("^.*__", "", .)))

Abund_strong_OTU_tax <- Abund_strong_2 %>%
  select(Genome, phylum)
sample_pop_all_list<- read.table("Sample_pop_meta/metadata_all_pop.txt", sep ="\t", header = TRUE)

dat_ggplot <- Abund_strong_2 %>%
  #select(-original_bin,-classification,-fastani_reference,fastani_reference_radius,-fastani_taxonomy,-fastani_ani,-fastani_af,-closest_placement_reference,-closest_placement_radius,-closest_placement_taxonomy,-closest_placement_ani,-closest_placement_af,-pplacer_taxonomy,-classification_method,-note,-other_related_references.genome_id.species_name.radius.ANI.AF.,-msa_percent,-translation_table,-red_value,-warnings,-completeness,-contamination,-length,-N50,-centrality, -fastani_reference_radius)
  #separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  select(Genome, starts_with(c("ERR","SRR"))) %>%
  pivot_longer(cols = -Genome, names_to = "sample_id", values_to = "counts" ) %>%
  left_join(Abund_strong_OTU_tax, by = "Genome") %>%
  left_join(sample_pop_all_list, by= "sample_id") %>%
  filter(!is.na(population)) %>%
  group_by(population,phylum) 


dat_ggplot$counts <- as.numeric(dat_ggplot$counts)
dat_ggplot$population <- factor(dat_ggplot$population, c("Hadza","Orang Asli","Nepali","Baka","Matses","Yanomami","Norman","Urban malaysian","California"))
dat_ggplot$phylum <- factor(dat_ggplot$phylum, c("Bacteroidota","Bacillota_A","Spirochaetota","Actinomycetota","Pseudomonadota","Verrucomicrobiota","Bacillota","Bacillota_B","Bacillota_C","Cyanobacteriota","Elusimicrobiota","Methanobacteriota","Desulfobacterota","Campylobacterota","Fibrobacterota","Myxococcota","Thermoplasmatota","Eremiobacterota","Planctomycetota","Fusobacteriota"))

write.csv(dat_ggplot, "data/abundance_phylum_pop_long.csv" ,row.names = FALSE)

# Installe le package si besoin
install.packages("scales")

# Charge le package
library(scales)


my_colors <- c(
  "#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#ffbb78",
  "#98df8a", "#aec7e8", "#ff9896", "#9467bd", "#c5b0d5",
  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#8B008B",
  "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
)


names(my_colors) <- c("Bacteroidota","Bacillota_A","Spirochaetota","Actinomycetota","Pseudomonadota","Verrucomicrobiota","Bacillota","Bacillota_B","Bacillota_C","Cyanobacteriota","Elusimicrobiota","Methanobacteriota","Desulfobacterota","Campylobacterota","Fibrobacterota","Myxococcota","Thermoplasmatota","Eremiobacterota","Planctomycetota","Fusobacteriota")

plot_AB<-ggplot(dat_ggplot, aes(x = sample_id, y=counts, fill= phylum ))+
  facet_grid(~ population, scales= "free_x", space = "free_x") +
  geom_bar(aes(fill = phylum), stat = "identity")+
  labs(y = "Relative abundance", x = "Sample") +
  theme( axis.text.x  =element_blank(), axis.ticks.x = element_blank(),strip.text.x = element_text(angle = 90, hjust = 0),
         strip.background = element_blank())+ 
  scale_fill_manual(values = my_colors)

plot_AB
ggsave(plot = plot_AB, "Figures/abundance_15mil_s95.png") 

dat_ggplot %>%
  group_by(sample_id) %>%
  summarize(total_abundance = sum(counts)) %>%
  print(n = Inf)  

#########alpha diveristy
library(vegan)
# Select colunnes
abund_cols <- Abund_strong_2 %>% select(starts_with(c("ERR","SRR"))) %>%
  mutate(across(starts_with(c("ERR","SRR")), as.numeric))

# Transposer to have the samples in the rows
abund_matrix <- abund_cols %>% 
  as.data.frame() %>%
  t() 
# verify is rows sums at 1
rowSums(abund_matrix)[1:5]

# Richness 
richness <- rowSums(abund_matrix > 0)

# Shannon 
shannon <- diversity(abund_matrix, index = "shannon")

# Simpson
simpson <- diversity(abund_matrix, index = "simpson")

alpha_div <- data.frame(
  sample_id = rownames(abund_matrix),
  richness = richness,
  shannon = shannon,
  simpson = simpson
)



head(alpha_div)
colnames(sample_pop_all_list)<- c("sample_id", "population")
alpha_div_pop <- merge(alpha_div, sample_pop_all_list, by ="sample_id")

shannon_alpha <- ggplot(alpha_div_pop, aes(x = population, y = shannon)) +
  geom_point() +
  theme_minimal() +
  labs(title = " Shannon diversity per population", y = "Shannon", x = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

simpson_alpha <- ggplot(alpha_div_pop, aes(x = population, y = simpson)) +
  geom_point() +
  theme_minimal() +
  labs(title = " Simpson diversity per population", y = "Simpson", x = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

richness_alpha <- ggplot(alpha_div_pop, aes(x = population, y = richness)) +
  geom_point() +
  theme_minimal() +
  labs(title = " Richness per population", y = "Richness", x = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = shannon_alpha, "Figures/shannon_alpha_15mil_s95.png") 
ggsave(plot = simpson_alpha, "Figures/simpson_alpha_15mil_s95.png") 
ggsave(plot = richness_alpha, "Figures/richness_alpha_15mil_s95.png") 

#########gamma diversity 
gamma_rich <- alpha_div_pop %>%
  group_by(population) %>%
  summarise(mean_richness = mean(richness, na.rm = TRUE))
gamma_shannon <- alpha_div_pop %>%
  group_by(population) %>%
  summarise(mean_shannon = mean(shannon, na.rm = TRUE))
gamma_simpson <- alpha_div_pop %>%
  group_by(population) %>%
  summarise(mean_simpson = mean(simpson, na.rm = TRUE))

gamma_all <- data.frame(
  Indice = c("Shannon", "Simpson", "Richness"),
  Valeur = c(gamma_shannon, gamma_shannon, gamma_rich)
)
gamma_all <- list(gamma_rich, gamma_shannon, gamma_simpson) %>% 
  reduce(full_join, by = "population")

# ggplot2 format (long)
gamma_all_long <- gamma_all %>%
  pivot_longer(cols = -population, 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("mean_richness", "mean_shannon", "mean_simpson"),
                         labels = c("Richness", "Shannon", "Simpson")))


# 2. Graphique combiné
# Création d'une palette de couleurs distinctes pour les populations
library(cowplot)
pop_colors <- brewer.pal(n = length(unique(gamma_all_long$population)), name = "Set3")


gamma_all_long$population<-factor(gamma_all_long$population, levels = c("Hadza","Baka","Matses", "Nepali", "Yanomami","California","Norman","Orang Asli", "Urban malaysian"))
richness_ggplot <- ggplot(gamma_all_long, aes(x= metric,y = value))+
  geom_point()
print(richness_ggplot)

#split dataset 

df_richness <- gamma_all_long %>% filter(metric == "Richness")
df_shannon <- gamma_all_long %>% filter(metric == "Shannon")
df_simpson <- gamma_all_long %>% filter(metric == "Simpson")

######beta diversity
abund_matrix_v2 <- na.omit(abund_matrix)
bray_dist <- vegdist(abund_matrix_v2, method = "bray")

####################Prevalence 

write.csv(Abund_strong_2, "./Table_phylum.csv")

# Load necessary library
library(dplyr)


# Calculate prevalence: fraction of samples where abundance > 0
Abund_strong_2$Prevalence <- rowSums(Abund_strong_2[, sampleid] > 0) / length(sampleid)

# Optional: View top OTUs by prevalence
Top_otu <-head(Abund_strong_2[order(-Abund_strong_2$Prevalence), c("Genome", "phylum", "species", "Prevalence")], n=20L)
write.csv(Abund_strong_2, "./Table_phylum.csv")


################subset aboundstrong 
df_phylum <- Abund_strong_2[,c("Genome","original_bin","phylum", "class","order","family", "genus", "species")]
write.csv(x = df_phylum, "./Name_species_phylum.csv")
write.csv(Abund_strong_2, "./Table_phylum_aboundance_prevalence.csv")



#####Stat 
# Required packages
library(dplyr)
library(tidyr)
library(rstatix)  # simplifie les tests pairwise


# Metrics to test
alpha_metrics <- c("richness", "shannon", "simpson")

# Empty list to store results
all_tests <- list()

alpha_div_pop_final <- alpha_div_pop %>%
  select( sample_id, richness, shannon, simpson, population)
# Loop over each metric and perform pairwise Wilcoxon tests
for (metric in alpha_metrics) {
  test_results <- alpha_div_pop_final %>%
    pairwise_wilcox_test(
      formula = as.formula(paste(metric, "~ population")),
      p.adjust.method = "fdr"
    ) %>%
   mutate(metric = metric)
  
  all_tests[[metric]] <- test_results
}


saveRDS(all_tests, "Results/P_alpha_adjusted_15mil_div.rds")




# Combine all results into a single dataframe
final_results <- bind_rows(all_tests)

# Show results
print(final_results)

## Figure2A and Figure S2
library(purrr)
library(ggsignif)
library(svglite)
install.packages("svglite")
mycols <- my_colors

select.metrics <- c("richness", "shannon", "simpson")
for(i in select.metrics){
  print(i)
  sumdata <- alpha_div_pop_final[, c("sample_id", "population", i)]
  colnames(sumdata)[3] <- "value"  # simplifie l'usage dans ggplot
  
 # urbpair <- list(c("Hadza", "California"),
       #         c("Hadza","Nepali"),
         #         c("Hadza","Baka"), 
          #        c("Baka", "California"))
  
  mytest <- sumdata %>%
    pairwise_wilcox_test(
      formula = value ~ population,
      p.adjust.method = "fdr"
    )
  
  # Extraire les paires significatives
  signif_pairs <- mytest %>% 
    filter(p.adj < 0.05) %>%
    mutate(pair = map2(group1, group2, c))  # crée une colonne de paires
  # Créer un dataframe pour la heatmap
  heatmap_df <- mytest %>%
    select(group1, group2, p.adj, p.adj.signif) %>%
    mutate(p.adj.signif = factor(p.adj.signif, 
                                 levels = c("ns", "*", "**", "***", "****")))  # ordre logique
  
  # # Position des étoiles : au-dessus du max + petits décalages
  # max_y <- max(sumdata$value, na.rm = TRUE)
  # signif_pairs <- signif_pairs %>%
  #   mutate(y.position = seq(0.6, 0.6 + 0.05 * (n() - 1), by = 0.05) * max_y)
  

  
 
  

  
  p <- ggplot(sumdata, aes(x = population, y = value)) + 
    theme_bw(base_size = 12, base_rect_size = 2, base_line_size = 0.8) +
    theme(panel.grid = element_blank()) +
    
    # Boxplots sans outliers
    geom_boxplot(size = 0.8, col = "grey50", outlier.shape = NA) +
    
    # # Étoiles de significativité : comparaisons entre populations
    # geom_signif(
    #   comparisons = signif_pairs$pair,
    #   annotations = signif_pairs$p.adj.signif,
    #   y_position = signif_pairs$y.position,
    #   size = 0.8,
    #   textsize = 4
    # ) +
    
    # Points individuels (jitter)
    geom_point(
      aes(colour = population, fill = population),
      position = position_jitter(width = 0.2, seed = 1),
      size = 2,
      shape = 21
    ) +
    
    # Thème et échelles
    scale_y_continuous(expand = c(0.07, 0)) +
    scale_color_manual(values = mycols) +
    scale_fill_manual(values = mycols) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank()) +
    ylab(i)
  
  
  p
  
  h <- ggplot(heatmap_df, aes(x = group1, y = group2, fill = p.adj.signif)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("ns" = "gray90", 
                 "*" = "#ffeda0", 
                 "**" = "#feb24c", 
                 "***" = "#f03b20", 
                 "****" = "#bd0026"),
      name = "Significance"
    ) +
    geom_text(aes(label = p.adj.signif), color = "black", size = 4) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(x = "", y = "", title = paste0("Pairwise significance (adjusted p-values)",i))
  h <- h +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5),  # marges générales
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.position = "bottom",                      # légende en bas
      legend.margin = margin(t = -5, unit = "pt")      # rapproche la légende de l’axe x
    ) +
    guides(fill = guide_legend(
      title.position = "top",
      title.hjust = 0,
      keywidth = 0.6,
      keyheight = 0.4
    )) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
  
  h
  
  
  # save
  if(i=="shannon"){
    filepath <- paste0("Figures/alpha_div_15mil_", i, ".svg")
    ggsave(filename = filepath, plot = p, device = "svg",width = 12, height = 12,units = "in", dpi="print", bg = "white")
    
    filepath <- paste0("Figures/alpha_div_15mil_", i, ".png")
    ggsave(filename = filepath, plot = p, device = "png",width = 12, height = 12,units = "in", dpi=600, bg = "white")
   
    filepath <- paste0("Figures/signif_alpha_div_15mil_", i, ".svg")
    ggsave(filename = filepath, plot = h, device = "svg",width = 8, height = 8,units = "in", dpi="print", bg = "white")
    
    filepath <- paste0("Figures/signif_alpha_div_15mil_", i, ".png")
    ggsave(filename = filepath, plot = h, device = "png",width = 8, height = 8,units = "in", dpi=600, bg = "white")
    # filepath <- paste0("Figures/alpha_div_15mil_", i, ".rds")
    # saveRDS(p, filepath)
  }else{
     filepath <- paste0("Figures/alpha_div_15mil_", i, ".svg")
     ggsave(filename = filepath, plot = p, device = "svg",width = 12, height = 12,units = "in", dpi="print", bg = "white")
    
    filepath <- paste0("Figures/alpha_div_15mil_", i, ".png")
    ggsave(filename = filepath, plot = p, device = "png",width = 12, height = 12,units = "in", dpi=600, bg = "white")
    
    filepath <- paste0("Figures/signif_alpha_div_15mil_", i, ".svg")
    ggsave(filename = filepath, plot = h, device = "svg",width = 8, height = 8,units = "in", dpi="print", bg = "white")
    
    filepath <- paste0("Figures/signif_alpha_div_15mil_", i, ".png")
    ggsave(filename = filepath, plot = h, device = "png",width = 8, height = 8,units = "in", dpi=600, bg = "white")
    # filepath <- paste0("Figures/alpha_div_15mil_", i, ".rds")
    # saveRDS(p, filepath)
  }
  
  print(p)
  print(h)
}

library(gridExtra)
lay <- rbind(c(rep(1,3),
               rep(2,3)))


# draw plot
diet.list <- grid.arrange(p,h,layout_matrix=lay)

lay <- rbind(c(rep(1,2), rep(2,2)))  # layout: 2 colonnes pour p, 2 pour h

grid.arrange(p, h,
             layout_matrix = lay,
             widths = c(1.5, 1.5, 1, 1))  # p = 5/7, h = 2/7


################# Kruskal Wallis 
# Load the necessary data
df <- dat_ggplot
df$counts <- as.numeric(df$counts)
# Calculate the total abundance by phylum, sample_id, and population
df_summary <- df %>%
  group_by(phylum, sample_id, population) %>%
  summarise(total_abundance = sum(counts)) %>%
  ungroup()



# Initialize lists to store results
kruskal.out <- list()
wilcox.out <- list()  # List to store Wilcoxon results
effect.out <- list()

# Perform Kruskal-Wallis tests for each phylum
for (i in unique(df_summary$phylum)) {  # Replace df_melted with df_summary if you already have a long-format dataframe
  tmp.df <- subset(df_summary, phylum == i)  # Filter by phylum
  
  # Kruskal-Wallis test (general test)
  krus.out <- kruskal.test(total_abundance ~ population, data = tmp.df)  
  
  # Calculate effect size with Kruskal-Wallis
  eff.out <- rstatix::kruskal_effsize(total_abundance ~ population, data = tmp.df, ci = TRUE, nboot = 100)  
  
  # Perform Wilcoxon pairwise test for specific comparisons between clusters
  pwc.out <- pairwise.wilcox.test(tmp.df$total_abundance, g = tmp.df$population, p.adjust.method = "BH")
  
  # Calculate effect size for Wilcoxon pairwise by passing the whole dataframe
  eff.out_wilcox <- wilcox_effsize(total_abundance ~ population, data = tmp.df)  # Correct use of data argument
  
  # Store results in the lists
  kruskal.out[[i]] <- krus.out
  effect.out[[i]] <- eff.out$effsize[1]  # Effect size from Kruskal-Wallis
  wilcox.out[[i]] <- list(pwc.out = pwc.out, effect_size = eff.out_wilcox$effsize)  # Add Wilcoxon effect size
}

# Create a summary table for Kruskal-Wallis results
krus.signif <- data.frame(
  phylum = names(kruskal.out),
  stat = unlist(lapply(kruskal.out, function(x) x$statistic)),
  pval = unlist(lapply(kruskal.out, function(x) x$p.value)),
  effsize = unlist(effect.out)
)

# Apply Benjamini-Hochberg correction for p-values
krus.signif$qval <- p.adjust(krus.signif$pval, method = "BH")

# Save Kruskal-Wallis results to a file
write.table(krus.signif, file = "output_files/Kruskal_Wallis_results.tsv", row.names = FALSE, sep = "\t")

# Filter and display significant phyla from Kruskal-Wallis results
krus.signif <- krus.signif[krus.signif$qval < 0.05, ]
krus.signif <- krus.signif[order(krus.signif$stat, decreasing = TRUE), ]
View(krus.signif)

# Wilcoxon test - Pairwise comparisons
# Initialize an empty data frame to store Wilcoxon results
wilcox.signif <- data.frame()

# Loop through each phylum in the Wilcoxon results
for (i in names(wilcox.out)) {
  out <- wilcox.out[[i]]$pwc.out$p.value
  effect_size <- wilcox.out[[i]]$effect_size  # Get effect size for Wilcoxon
  
  # Convert the p-value matrix into a long format using melt
  wilcox.melted <- melt(out)  # This will melt the p-values into long format
  
  # Align the effect size correctly with the melted data
  # Repeat the effect size for each comparison between clusters (number of rows in wilcox.melted)
  wilcox.melted$effect_size <- rep(effect_size, length.out = nrow(wilcox.melted))  # Correctly repeat effect size
  
  # Add the 'family' column (phylum) to the melted data
  wilcox.melted$family <- i
  
  # Combine the results into the final dataframe
  wilcox.signif <- rbind(wilcox.signif, wilcox.melted)
}

# Remove NA values (optional, but ensures we don't have incomplete data)
wilcox.signif <- na.omit(wilcox.signif)

# Filter significant results with p-value < 0.05
#wilcox.signif <- wilcox.signif[wilcox.signif$value < 0.05, ]

# Sort the results by the significance level (p-value)
#wilcox.signif <- wilcox.signif[order(wilcox.signif$value), ]

# Save the Wilcoxon results to a file
write.table(wilcox.signif, file = "output_files/Wilcoxon_Family_vs_Cluster_Pairwise.tsv", row.names = FALSE, sep = "\t")

wilcox.signif_test <- wilcox.signif %>%
  mutate(Var1 = recode(Var1, 
                             "Urban malaysian" = "UB", 
                             "Orang Asli" = "OA")) %>%
  mutate(Var2 = recode(Var2, 
                       "Urban malaysian" = "UB", 
                       "Orang Asli" = "OA"))
#####plot 
#####plot 
a <- ggplot(wilcox.signif_test, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color="white", linewidth=0.1) +  # Utiliser geom_tile pour créer la heatmap
  scale_fill_viridis( breaks = c(0, 0.05, 1),  # Définir les breaks pour la légende
                      labels = c("0","0.05", "1")) + # Étiquettes de la légende) +  # Couleurs de la heatmap
  facet_wrap(~family) +  # Créer une heatmap pour chaque phylum
  theme_minimal() +  # Thème minimal
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotation des labels des axes x
  labs(title = "Heatmap of Pairwise Comparisons by Phylum", 
       x = "Population ", 
       y = "Population", 
       fill = "P-value")  + theme_tufte(base_family="Helvetica") +
  theme(legend.title=element_text(size=8)) + theme(panel.spacing.y=unit(0.1, "cm")) +
  theme(panel.margin.x=unit(0.1, "cm"))+ 
  theme(legend.key.width=unit(0.5, "cm"))+
   theme(legend.key.size=unit(1, "cm")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 15)))
a
ggsave("Figures/heatmap_pairwise_phylum_pop.pdf", plot = a, device = "pdf", width = 10, height = 8)


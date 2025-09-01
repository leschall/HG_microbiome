library(tidyverse)
library(ape)
library(vegan)
library(geosphere)
library(maps)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
install.packages("geosphere")
library(geosphere)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
library("ggtree")
# 
# c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
#          "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
# read_tsv("../Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(! study_name == "ThomasAM_2019_c" ) -> Data
# unique(Data$Country) -> C_name
# c25[1:length(C_name)] -> C_colors
# names(C_colors) = C_name

library(dplyr)
library(stringr)

# Exemple de dataframe
def_meta <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ";")
mantel_df <- def_meta

# Traitement
df_clean <- mantel_df %>%
  mutate(
    # Supprimer les espaces multiples
    coord = str_squish(Lat.Long),
    # Extraire latitude, longitude, et leurs directions
    latitude = as.numeric(str_extract(coord, "^[0-9\\.]+")),
    lat_dir = str_extract(coord, "(?<= )[NS](?= )"),
    longitude = as.numeric(str_extract(coord, "(?<= )[0-9\\.]+(?= [EW])")),
    lon_dir = str_extract(coord, "[EW]$")
  ) %>%
  # Appliquer les signes négatifs si direction est S ou W
  mutate(
    latitude = if_else(lat_dir == "S", -latitude, latitude),
    longitude = if_else(lon_dir == "W", -longitude, longitude)
  ) %>%
  select(-lat_dir, -lon_dir)

print(df_clean)

#distance matrix 
#longitude and latitude 
geo = data.frame(df_clean$longitude, df_clean$latitude)
rownames(geo) <- rownames(df_clean)  # ou noms des samples


# 2. Lire ton arbre (Newick)
tree <- read.tree(file ="data/Phylogeny_geography/1002_1_mugsy_tree.treefile")
tree_test <- read.tree(file = "data/Phylogeny_geography/1002_1_parsnp_tree.treefile")
ggtree(tree_test) + geom_tiplab()
ggtree(tree) + geom_tiplab()


tree[["tip.label"]] <- gsub("_bin$", "", tree[["tip.label"]])

#Merge data and tree information
ID_anal <- tree[["tip.label"]]
geo_anal <- geo[rownames(geo) %in% ID_anal, ]
#keep.tip(Tree, Data_anal$ID_anal) -> Tree
order_indices <- match(tree$tip.label, as.character(rownames(geo_anal)))
geo_anal[order_indices,] -> geo_anal
all(rownames(geo_anal) == ID_anal)


#geographic data frame - haversine distance 
d.geo = distm(geo_anal, fun = distHaversine)
rownames(d.geo) <- rownames(geo_anal)
colnames(d.geo) <- rownames(geo_anal)
dist.geo = as.dist(d.geo)
#tree distance 
ggtree(tree) + geom_tiplab()
dist_tree <- cophenetic.phylo(tree)
dist_tree <- as.dist(dist_tree)



write.table(dist.geo, "Results/Phylogeny_geography_Mandel/dist.geo.txt")
write.table(dist_tree, "Results/Phylogeny_geography_Mandel/dist.tree.txt")

##mantel tesdist_tree##mantel test 
Perm =9999
mantel( dist_tree, dist.geo, method = "pearson", permutations=Perm, parallel=1 ) -> Test
tibble(Rho= Test$statistic, P=Test$signif, Permutations= Test$permutations) -> Results


#############################################################
####metadata





# Chargement des packages
library(tidyverse)
library(ape)
library(vegan)
library(geosphere)

# Chemins
tree_dir <- "data/Phylogeny_geography/Alltrees"
out_dir <- "Results/Phylogeny_geography_Mandel"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Lecture des métadonnées
mantel_df <- input.metadata

# Nettoyage des coordonnées
df_clean <- mantel_df %>%
  mutate(coord = str_squish(Lat.Long),
         latitude = as.numeric(str_extract(coord, "^[0-9\\.]+")),
         lat_dir = str_extract(coord, "(?<= )[NS](?= )"),
         longitude = as.numeric(str_extract(coord, "(?<= )[0-9\\.]+(?= [EW])")),
         lon_dir = str_extract(coord, "[EW]$"),
         latitude = if_else(lat_dir == "S", -latitude, latitude),
         longitude = if_else(lon_dir == "W", -longitude, longitude)) %>%
  select(sample_id, latitude, longitude)

# Mise en forme de la matrice géographique
geo <- df_clean %>% select(longitude, latitude)
rownames(geo) <- df_clean$sample_id

# Fonction pour extraire le préfixe identifiant de l'arbre (ex: "439_1")
SGB_from_tree <- function(filename){
  split_filename <- unlist(strsplit(basename(filename), "_"))
  return(paste(split_filename[1], split_filename[2], sep = "_"))
}
SGB_from_tree("data/Phylogeny_geography/Alltrees/936_45_mugsy_tree.treefile")

# Fonction qui fait l’analyse complète pour un arbre
Run_analysis <- function(tree_file, Perm = 9999) {
  tree <- read.tree(tree_file)
  tree$tip.label <- gsub("_bin$", "", tree$tip.label)
  ID_anal <- tree$tip.label
  
  geo_anal <- geo[rownames(geo) %in% ID_anal, ]
  if (nrow(geo_anal) < 3) return(NULL)
  
  geo_anal <- geo_anal[match(ID_anal, rownames(geo_anal)), , drop = FALSE]
  if (any(is.na(rownames(geo_anal)))) {
    cat("Warning: unmatched samples for", basename(tree_file), "\n")
    return(NULL)
  }
  
  d.geo <- distm(geo_anal, fun = distHaversine)
  rownames(d.geo) <- rownames(geo_anal)
  colnames(d.geo) <- rownames(geo_anal)
  dist.geo <- as.dist(d.geo)
  
  dist.tree <- as.dist(cophenetic.phylo(tree))
  
  Test <- mantel(dist.tree, dist.geo, method = "pearson", permutations = Perm)
  tibble(Rho = Test$statistic, P = Test$signif, Permutations = Test$permutations)
}

# Liste des fichiers d’arbres
tree_files <- list.files(tree_dir, pattern = "\\.treefile$", full.names = TRUE)

# Résultat global
Results <- tibble()
Perm= 9999
# Boucle sur tous les arbres
for (tree_file in tree_files) {
  SGB <- SGB_from_tree(tree_file)
  Out <- file.path(out_dir, paste0(SGB, ".tsv"))
  tree <- read.tree(tree_file)
  tree$tip.label <- gsub("_bin$", "", tree$tip.label)
  n_tips <- length(tree$tip.label)
  
  if (file.exists(Out)) {
    cat(SGB, ": Already done\n")
    next
  }
  
  cat("Processing", SGB, "\n")
  result <- Run_analysis(tree_file)
  
  if (is.null(result)) {
    cat(SGB, ": Skipped (bad data)\n")
    next
  }
  
  result <- result %>% mutate(SGB = SGB, N_samples = n_tips)
  write_tsv(result, Out)
  Results <- bind_rows(Results, result)
}

# Sauvegarde finale
write_tsv(Results, file.path(out_dir, "All_cor_2.tsv"))

## correction 
adj_pval <- read.table(file = "Results/Phylogeny_geography_Mandel/All_cor.tsv", header=T)

adj_pval$P_adj <- p.adjust(p = adj_pval$P_rho, method = "fdr")  # équivalent à method = "BH"
write_csv(adj_pval, "Results/Phylogeny_geography_Mandel/All_cor_adj.csv")

###import metadata plot 




cluster_classification<-read.csv( "data/cluster_strains_final.csv ", header = T)
cluster_classification <- cluster_classification %>% 
  select(secondary_cluster, phylum, class, order, family, genus, species) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>% 
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  distinct(secondary_cluster, .keep_all = TRUE)
cluster_classification <- cluster_classification %>% rename(cluster = secondary_cluster)
merge_test<- read.csv(file = "Results/Phylogeny_geography_Mandel/All_cor_adj.csv")

df_plot_taxon <- merge(merge_test, cluster_classification, by = "cluster")  # full outer join
write.csv(df_plot_taxon, "Results/Phylogeny_geography_Mandel/roh_mantel_taxon.csv", row.names = F)
library(ggplot2)

ggplot(Results, aes(x = Rho)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  labs(title = "Distribution des coefficients de Mantel (Rho)",
       x = "Coefficient de Mantel (Rho)", y = "Nombre d'arbres") +
  theme_minimal()


ggplot(Results, aes(x = Rho, y = -log10(P))) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Corrélation géographique vs significativité",
       x = "Rho (Mantel)", y = "-log10(p-value)") +
  theme_minimal() +
  ylim(-0.5 , 4.5) 

ggplot(df_plot_taxon, aes(x = family, y = Rho, fill = family)) +
  geom_boxplot() +
  labs(title = "Comparaison des corrélations par taxon", y = "Rho (Mantel)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####PLot
library(tidytext) 
mantel_res <-read.csv(file = "Results/Phylogeny_geography_Mandel/roh_mantel_taxon.csv", sep = ";")
mantel_res$phylum <- factor(mantel_res$phylum, c("Bacteroidota","Bacillota_A","Spirochaetota","Actinomycetota","Pseudomonadota","Verrucomicrobiota","Bacillota","Bacillota_B","Bacillota_C","Cyanobacteriota","Elusimicrobiota","Methanobacteriota","Desulfobacterota","Campylobacterota","Fibrobacterota","Myxococcota","Thermoplasmatota","Eremiobacterota","Planctomycetota","Fusobacteriota"))
# order data

mantel_res <- mantel_res[order(mantel_res$Rho,decreasing = F),]
mantel_res <- mantel_res %>%
  mutate(species_ord = reorder_within(species, Rho, phylum))
a <- ggplot(mantel_res, aes(x= Rho,y= species_ord))+ 
  theme_minimal() + 
  theme(plot.margin = margin(0,0,21,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid=element_blank(),
        panel.grid.major.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_rect(fill="white",linewidth = 0),
        text=element_text(size=12), axis.text=element_text(size=12),
        #axis.text.y = element_text(angle=0, hjust = 1, size=8),
        axis.text.x=element_text(colour = "black", size=8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y= element_blank())+ 
  scale_color_manual(values = "pink") +
  scale_x_continuous(breaks = seq(-0.1,1, by=0.05),limits = c(-0.1,1), expand = c(0.02,0.06))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Rho") +
  labs(title="Significant Rho (pval<0.05)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  
  geom_vline(xintercept = c(0, 0.3), colour = "red", linetype = "dashed", linewidth = 0.7, alpha = 0.9) +
  geom_point(size=3, colour="pink") +
  geom_text(aes(label = species), nudge_x = 0.1, size=3, col="grey40", hjust=0, vjust=0.2) +
  geom_point(aes(col= "pink"),size=3) +
  geom_text(aes(label=species), angle=0, nudge_x = 0.1, 
            size=3, col="grey40", hjust=0, vjust=0.2)
a
ggsave("Results/Phylogeny_geography_Mandel/rho_representation.png",device = "png", a, width = 10, height = 20, dpi = 300)



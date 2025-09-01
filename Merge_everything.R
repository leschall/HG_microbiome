####Merge files 

Summary_gini_results <- read.csv(file = "Results/Rtapas/summary_gini_results.csv", header= TRUE, sep =",")
paco_summary_results <- read.csv(file = "Results/Paco/paco_summary_results.csv", header= TRUE, sep =",")
Anpan_results <- read.csv(file = "Results/Anpan/anpan_effects_reshaped.csv", header= TRUE, sep =",")
phylogeogprahy_results <-  read.csv(file = "Results/Phylogeny_geography_Mandel/All_cor_adj.csv", header= TRUE)

Anpan_results <- Anpan_results %>% 
  rename(cluster = cluster_id)

all_results <- phylogeogprahy_results %>%
  left_join(Summary_gini_results, by="cluster")

all_results <- all_results %>% 
  left_join(paco_summary_results, by = "cluster")

all_results <- all_results %>%
  left_join(Anpan_results, by = "cluster")

###Phylogeny 
cluster_classification<-read.csv( "data/cluster_strains_final.csv ", header = T)

cluster_classification <- cluster_classification %>% 
  select(secondary_cluster, phylum, class, order, family, genus, species) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>% 
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  distinct(secondary_cluster, .keep_all = TRUE)
cluster_classification <- cluster_classification %>% rename(cluster = secondary_cluster)
write.csv(cluster_classification, "data/cluster_classification_final.csv" )

all_results <- all_results %>%
  left_join(cluster_classification, by= "cluster")

write.csv( all_results, "Results/all_results_anpan_rtapas_paco_phylo.csv", row.names = FALSE)



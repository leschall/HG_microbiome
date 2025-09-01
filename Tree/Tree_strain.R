###########Tree strains 

library("dplyr")
library("tidyr")
cluster_tree <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Tree/Cdb.csv")
Human_fecal_metadata <-read.table("/Users/lucie/Documents/Cambridge/Analysis/Tree/HumanFecal_modernHG_draft1_metadata.tsv", header= T)
sample_pop_all <- read.csv(file = "Sample_pop_meta/sample_pop_all.list.csv", header=T, sep = ";")
cluster_one <- cluster_tree %>%
  group_by(primary_cluster) %>%
  count(primary_cluster) %>%
  filter(n > 10) 
cluster_second <- cluster_tree %>%
  filter(cluster_tree$primary_cluster %in% cluster_one$primary_cluster) %>%
  group_by(secondary_cluster) %>%
  count(secondary_cluster) %>%
  filter(n > 10)

cluster_df <- cluster_tree %>%
  separate(genome, into = c("sampleid", "bin"), sep = "_") %>%
  left_join(sample_pop_all, by = "sampleid") %>%
  filter(secondary_cluster %in% cluster_second$secondary_cluster)%>%
  group_by(secondary_cluster) %>%
  filter(n_distinct(population) >= 3) %>%
  ungroup()%>%
  unite("original_bin", sampleid, bin, sep = "_") %>%
  mutate(original_bin = sub("\\.fa$", "", original_bin)) %>% 
  left_join(Human_fecal_metadata, by="original_bin")

cluster_species <- cluster_df %>%
  mutate(classification = na_if(classification, "NA "))%>%
  group_by(secondary_cluster) %>%
  mutate(classification = first(na.omit(classification))) %>%
  ungroup() %>% 
  select(c("original_bin","threshold","cluster_method","comparison_algorithm","primary_cluster",
           "original_secondary_cluster","secondary_cluster","population","Age_bin","Age","Gender",
           "Country","Lat.Long","Village","SampleName","Lifestyle","State","BMI","Author","Year",
           "classification","completeness","contamination","length","N50")) %>%
  separate(classification, into = c("domain", "phylum", "class","order","family","genus","species"), sep = ";") %>%
  mutate(genus = sub("^g__", "", genus))


write.csv(cluster_species, "data/cluster_strains_final.csv ", row.names = FALSE)
write.table(cluster_df, "data/cluster_strains.txt", sep ="\t", row.names = FALSE)

reads_count <- read.csv("data/N_count_reads_clean_sample.csv", sep=";")
metadata_all <- merge(metadata, reads_count, by ="sample_id")
write.csv(metadata_all, "Sample_pop_meta/metadata_pop_all_final.csv", row.names = FALSE )

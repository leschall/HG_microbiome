# -------------------------------
# Libraries
# -------------------------------
library(tidyverse)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")

library(fgsea)

# -------------------------------
# Load data
# -------------------------------
df <- read_csv("Results/Phylogeny_geography_Mandel/roh_mantel_taxon.csv", sep=";")
df <-read.csv(file = "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/significant_species.csv" , sep =";")
df <- read_csv("Results/Anpan/anpan_effects_reshaped.csv")
taxonomy <- read_csv("data/cluster_classification_final.csv")
df <- df %>% 
  rename(cluster = cluster_id) %>%
  left_join(taxonomy, by ="cluster")

# -------------------------------
# Define geographic association criteria
# -------------------------------
df <- df %>%
  mutate(Geo_associated = (Rho > 0.3 & P_rho < 0.05))
# -------------------------------
# Define anpan association criteria
# -------------------------------
df <- df %>%
  mutate(phylo_associated = (base_fit_elpd_diff < 0 & (base_fit_elpd_diff + 2* base_fit_se_diff)< 0))

# -------------------------------
# Define maaslin association criteria
# -------------------------------
df <- df %>%
  mutate(specie_associated = ( coef > 0 ))

# -------------------------------
# FUNCTION: Fisher enrichment
# -------------------------------
run_fisher_enrichment <- function(df, tax_level) {
  results <- tibble()
  subset_data <- df %>% filter(specie_associated)
  
  for (taxon in unique(df[[tax_level]])) {
    total_in_taxon <- sum(df[[tax_level]] == taxon)
    
    # Skip taxons with fewer than 5 SGBs
    if (total_in_taxon < 5) next
    
    sig_in_taxon <- sum(subset_data[[tax_level]] == taxon)
    nonsig_in_taxon <- total_in_taxon - sig_in_taxon
    
    sig_not_taxon <- sum(subset_data[[tax_level]] != taxon)
    total_not_taxon <- sum(df[[tax_level]] != taxon)
    nonsig_not_taxon <- total_not_taxon - sig_not_taxon
    
    contingency <- matrix(c(sig_in_taxon, nonsig_in_taxon,
                            sig_not_taxon, nonsig_not_taxon),
                          nrow = 2)
    
    test <- fisher.test(contingency, alternative = "greater")
    
    results <- bind_rows(results, tibble(
      Taxon = taxon,
      P = test$p.value,
      N_taxa = total_in_taxon,
      Sig_taxa = sig_in_taxon
    ))
  }
  
  results <- results %>%
    mutate(
      FDR = p.adjust(P, method = "fdr"),
      Taxonomic_level = tax_level
    )
  
  return(results)
}


# -------------------------------
# Run Fisher tests on all levels
# -------------------------------
fisher_all <- bind_rows(
  run_fisher_enrichment(df, "phylum"),
  run_fisher_enrichment(df, "class"),
  run_fisher_enrichment(df, "order"),
  run_fisher_enrichment(df, "family"),
  run_fisher_enrichment(df, "genus")
)

write_csv(fisher_all, "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop//Results_enrichment_Fisher.csv")

# -------------------------------
# RANK-BASED ENRICHMENT (fgsea)
# -------------------------------


if (is.character(df$coef)) {
  df <- df %>% mutate(coef = as.numeric(gsub(",", ".", coef)))
}
df_rho <- df %>% arrange(desc(coef))
ranks <- df_rho$coef
names(ranks) <- df_rho$feature 


#########ANNOYING DF 
df_clean <- df %>%
  mutate(feature = trimws(feature)) %>%
  filter(!is.na(feature), !is.na(coef), is.finite(coef)) %>%
  group_by(feature) %>%
  slice_max(order_by = abs(coef), n = 1, with_ties = FALSE) %>%
  ungroup()

family_sets <- df_clean %>%
  filter(!is.na(family), str_detect(family, "aceae")) %>%
  group_by(family) %>%
  summarise(features = list(unique(feature)), N = n(), .groups = "drop") %>%
  filter(N >= 5)

genus_sets <- df_clean %>%
  filter(!is.na(genus)) %>%
  group_by(genus) %>%
  summarise(features = list(unique(feature)), N = n(), .groups = "drop") %>%
  filter(N >= 5)

# restrict to the features present in ranks
genus_sets$features <- lapply(genus_sets$features, function(v) intersect(v, names(ranks)))
genus_sets <- genus_sets[lengths(genus_sets$features) > 0, ]
taxa_list <- setNames(genus_sets$features, genus_sets$genus)

ranks <- setNames(df_clean$coef, df_clean$feature)
fgseaRes <- fgsea(pathways = taxa_list, stats = ranks) %>% arrange(padj)


################





# Only families (ends with 'aceae') and not NA
family_sets <- df %>%
  filter(!is.na(family), str_detect(family, "aceae")) %>%
  group_by(family) %>%
  summarise(features = list(feature), N = n()) %>%
  filter(N >= 5)

#Only phylum  and not NA
# phylum_sets <- df %>%
#   filter(!is.na(phylum)) %>%
#   group_by(phylum) %>%
#   summarise(clusters = list(cluster), N = n()) %>%
#   filter(N >= 5)

# Convert to named list
taxa_list <- setNames(family_sets$features, family_sets$family)

# -------------------------------------
# Run fgsea only on selected families
# -------------------------------------
fgseaRes <- fgsea(pathways = taxa_list, stats = ranks)
fgseaRes <- fgseaRes %>% arrange(padj)

write_csv(fgseaRes, "Results/Phylogeny_geography_Mandel/Results_enrichment_fgsea.csv")

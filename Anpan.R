####ANPAN

install.packages(c("ape", 
                   "data.table",
                   "dplyr", 
                   "fastglm",
                   "furrr", 
                   "ggdendro",
                   "ggnewscale",
                   "ggplot2",
                   "loo",
                   "patchwork",
                   "phylogram",
                   "posterior",
                   "progressr",
                   "purrr", 
                   "R.utils",
                   "remotes",
                   "stringr",
                   "tibble",
                   "tidyselect")) # add Ncpus = 4 to go faster

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
remotes::install_github("biobakery/anpan")
library(anpan)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ape)
library(stringr)




df <- read.delim("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/data/Table_final_OTU.txt", header = TRUE, sep= "\t")
meta_rtapas <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ";")

cluster_tree <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Tree/Cdb.csv")
tree_strain_test <- read.tree(file = "data/Phylogeny_geography/Alltrees/816_24_mugsy_tree.treefile")

### filter info meta_anapan 
# separate by cluster 
cluster_anpan<-read.csv( "data/cluster_strains_final.csv ", header = T)
cluster_anpan <- merge(cluster_anpan, meta_rtapas, by = "original_bin")
cluster_anpan <- cluster_anpan %>% 
  select(original_bin,secondary_cluster,population,Age_bin,Gender, Country, Lifestyle, Author ) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>% 
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  mutate(sample_id = original_bin)%>%
  mutate(sample_id= str_remove(sample_id, "\\..*$"))

df_subset <- def_meta %>% select(sample_id, Geographic_region, N_reads_clean) %>% 
  rename(sample_ide = sample_id)
  
cluster_anpan <- cluster_anpan %>%
  mutate(sample_ide = str_remove(original_bin, "_bin\\..*$")) %>%
  left_join(df_subset, by = "sample_ide")
  
#Filter les datas
anpan_filtré <- cluster_anpan %>%
  filter(secondary_cluster == "816_24")
anpan_filtré$Lifestyle = factor(anpan_filtré$Lifestyle, levels = c("Industrial", "H-G"))
anpan_filtré$Geographic_region <- factor(anpan_filtré$Geographic_region, levels =c("America","Africa","Asia"))
######tree

plot_outcome_tree(tree_file = tree_strain_test, meta_file = anpan_filtré, covariates=NULL, outcome="Lifestyle")

####Run anpan 

result_test <- anpan_pglmm(
  meta_file = anpan_filtré ,
  tree_file = tree_strain_test,
  outcome = "Lifestyle",
  covariates = "Geographic_region",
  out_dir = "Results/Anpan",
  family = "binomial",
  save_object = FALSE,
  parallel_chains = 4,
  reg_noise = FALSE
)




draws <- result$pglmm_fit$draws()
oula <- summarise_draws(draws)
write.csv(oula, "Results/Anpan/oula.csv")


##################################

#### ANPAN FULL PIPELINE WITH LOOP + PLOT SAVING

# Install packages
install.packages(c("ape", "data.table", "dplyr", "fastglm", "furrr", "ggdendro", "ggnewscale",
                   "ggplot2", "loo", "patchwork", "phylogram", "posterior", "progressr",
                   "purrr", "R.utils", "remotes", "stringr", "tibble", "tidyselect"))

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2)

# Load libraries
remotes::install_github("biobakery/anpan")
library(anpan)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ape)
library(stringr)
library(posterior)

# Load input data
df <- read.delim("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/data/Table_final_OTU.txt", sep = "\t")
meta_rtapas <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", sep = ";")
cluster_tree <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Tree/Cdb.csv")
cluster_anpan <- read.csv("data/cluster_strains_final.csv ", header = TRUE)

# Merge and clean metadata
cluster_anpan <- cluster_anpan %>%
  mutate(sample_ide = str_remove(original_bin, "_bin\\..*$"))
# Add additional metadata
df_subset <- meta_rtapas %>%
  select(sample_id, Geographic_region, N_reads_clean) %>%
  rename(sample_ide = sample_id)

cluster_anpan <- cluster_anpan %>%
  left_join(df_subset, by = "sample_ide")
cluster_anpan<- cluster_anpan %>%
  select(original_bin, secondary_cluster, population, Age_bin, Gender, Country, Lifestyle, Author, Geographic_region, N_reads_clean, sample_ide) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>%
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  mutate(sample_id = str_remove(original_bin, "\\..*$"))


# Define key variables
tree_dir <- "data/Phylogeny_geography/Alltrees"
tree_files <- list.files(tree_dir, pattern = "*.treefile", full.names = TRUE)
outcome_var <- "Lifestyle"
covariate_var <- c("Geographic_region")

# Set factor levels
cluster_anpan$Lifestyle <- factor(cluster_anpan$Lifestyle, levels = c("Industrial", "H-G"))
pop.ord <- as.character(unique(cluster_anpan$population))
n_basap <- which(pop.ord=="California")
pop.ord <- c(pop.ord[n_basap],pop.ord[-n_basap])
cluster_anpan$population <- factor(cluster_anpan$population, pop.ord)
cluster_anpan$Age_bin <- factor(cluster_anpan$Age_bin, levels = c("18 and above", "14 and below"))
cluster_anpan$Country <- factor(cluster_anpan$Country, levels =c("USA","Tanzania","Nepal","Brazil","Cameroon","Peru","Malaysia"))
cluster_anpan$Author <- factor(cluster_anpan$Author, levels =c("Carter","Conteville","Rampelli","Obregon-Tito","Tee"))
cluster_anpan$Geographic_region <- factor(cluster_anpan$Geographic_region, levels =c("America","Africa","Asia"))

# Initialize empty dataframe to store all beta and sigma_phylo rows
all_effects <- data.frame()
all_loo <- data.frame()
# Loop over each tree file
for (tree_path in tree_files) {
  try({
    # Extract cluster ID from filename
    cluster_id <- str_extract(basename(tree_path), "^[^_]+_[^_]+")
    
    # Load the tree
    tree <- read.tree(tree_path)
    
    # Filter metadata for current cluster
    anpan_filtré <- cluster_anpan %>% filter(secondary_cluster == cluster_id)
    if (nrow(anpan_filtré) < 2) {
      message(paste("Skipping", cluster_id, "- not enough samples"))
      next
    }
    
    # Create output folders
    out_dir <- file.path("Results/Anpan/Fixed_random/Lifestyle_count__geo_ran_author_pop", cluster_id)
    plot_dir <- file.path(out_dir, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Run ANPAN
    result <- anpan_pglmm(
      meta_file = anpan_filtré,
      tree_file = tree,
      outcome = outcome_var,
      covariates = covariate_var,
      out_dir = out_dir,
      family = "binomial",
      save_object = FALSE,
      parallel_chains = 4,
      reg_noise = FALSE
    )
    
    # Save main plot
    ggsave(
      filename = file.path(plot_dir, paste0("AnpanPlot_", outcome_var, "_", cluster_id, ".png")),
      plot = result$plot,
      width = 10, height = 8, dpi = 300
    )
    
    # Save summary of posterior
    draws <- result$pglmm_fit$draws()
    oula <- summarise_draws(draws)
    
    csv_name <- paste0("Oula_", outcome_var, "_", covariate_var, ".csv")
    write.csv(oula, file.path(plot_dir, csv_name), row.names = FALSE)
    
    # Filter effects of interest and store with cluster ID
    selected_effects <- oula %>%
      filter(str_detect(variable, "beta") | str_detect(variable, "sigma_phylo") | str_detect(variable, "centered_cov_intercept")) %>%
      mutate(cluster_id = cluster_id)
    
    all_effects <- bind_rows(all_effects, selected_effects)
    
    # Force into a data frame (if it is not already one)
    loo_comparison <- as.data.frame(result[["loo"]][["comparison"]])
    
    # dd model names (currently in the rownames) as a ‘model’ column
    loo_comparison$model <- rownames(loo_comparison)
    loo_comparison <- loo_comparison %>%
      relocate(model, .before = everything())
    loo_comparison <- loo_comparison %>%
      mutate(cluster_id = cluster_id)
    
    all_loo <- bind_rows(all_loo, loo_comparison)
    
    
    message(paste("Finished cluster:", cluster_id))
    
  }, silent = FALSE)
}

write.csv(all_effects, "Results/Anpan/anpan_all_effects.csv", row.names = FALSE)
write.csv(all_loo, "Results/Anpan/anpan_all_loo.csv", row.names = FALSE)
#####Transform to a large format 
library(dplyr)     # For data manipulation
library(tidyr)     # For reshaping data (pivoting)
library(readr)     # For reading CSV files

# Convert to long format 
all_effects_long <- all_effects %>%
  pivot_longer(cols = c(mean, median, sd, mad, q5, q95, rhat, ess_bulk, ess_tail),
               names_to = "statistic",
               values_to = "value")
all_effects_long <- all_effects_long %>%
  mutate(variable_stat = paste0(variable, "_", statistic))

all_effects_wide <- all_effects_long %>%
  select(cluster_id, variable_stat, value) %>%
  pivot_wider(names_from = variable_stat, values_from = value)

all_loo_long <- all_loo %>%
  select(model,elpd_diff,se_diff, cluster_id) 

write.csv(all_loo_long, "Results/Anpan/loo.csv")


# Pivot long → large pour avoir une colonne par combinaison model × variable
loo_wide <- all_loo_long %>%
  select(cluster_id, model, elpd_diff, se_diff) %>%
  pivot_wider(
    names_from = model,
    values_from = c(elpd_diff, se_diff),
    names_sep = "_"
  ) %>%
  # Rename the columns to make them more readable
  rename(
    base_fit_elpd_diff = elpd_diff_base_fit,
    pglmm_fit_elpd_diff = elpd_diff_pglmm_fit,
    base_fit_se_diff = se_diff_base_fit,
    pglmm_fit_se_diff = se_diff_pglmm_fit
  )


final_df <- left_join(all_effects_wide, loo_wide, by = "cluster_id") #%>%
  #filter(model != "pglmm_fit")

write_csv(final_df, "Results/Anpan/anpan_effects_reshaped.csv")

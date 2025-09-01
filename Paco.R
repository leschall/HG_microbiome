library(ape)
library(dplyr)
library(tidyr)
library(Rtapas)
library(stringr)
library(paco)

# Load metadata and input files
df <- read.delim("data/Table_final_OTU.txt", header = TRUE, sep= "\t")
meta_rtapas <- read.csv("Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ";")


# Load and process the host tree (human populations)
tree_human <- "(Baka,(Hadza,((California,Norman),(OrangAsli,(UrbanMalaysians,(Nepali,(Matses, Yanomami)))))));"
tree <- ladderize(compute.brlen(read.tree(text = tree_human), method = 1), right = FALSE)

# Load cluster assignments and clean IDs
cluster_ratapas <- read.csv("data/cluster_strains_final.csv ", header = TRUE)
cluster_rtapas <- cluster_ratapas %>% 
  select(original_bin, secondary_cluster) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>% 
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  mutate(sample_id = str_remove(original_bin, "_.*$"))

# Prepare relative abundance table (CPM-normalized)
reads <- df[,sampleid]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(comm.df = reads)
dim(reads)
input.rtapas <- as.data.frame(t(reads))
#input.rtapas <- input.comm[sampleid,]
input.rtapas$sample_id <- rownames(input.rtapas)


# Format metadata
input_meta_rtapa <- meta_rtapas %>%
  select(sample_id, population) %>%
  mutate(population = case_when(
    population == "Urban malaysian" ~ "UrbanMalaysians",
    population == "Orang Asli" ~ "OrangAsli",
    TRUE ~ population
  ))

# Merge metadata with abundance data
merged_rtapa <- merge(input.rtapas, input_meta_rtapa, by = "sample_id")


# Initialize global results dataframe
global_results <- data.frame()

# Loop over each microbial cluster
for (cluster_id in unique(cluster_rtapas$secondary_cluster)) {
  message("Processing cluster: ", cluster_id)
  
  cluster_sub <- cluster_rtapas %>% filter(secondary_cluster == cluster_id)
  sampleid <- cluster_sub$original_bin
  
  if (length(sampleid) < 5) next  # Skip small clusters
  
  strain_tree_file <- paste0("data/Phylogeny_geography/Alltrees/", cluster_id, "_mugsy_tree.treefile")
  if (!file.exists(strain_tree_file)) next
  tree_strain <- read.tree(file = strain_tree_file)
  
  # Create host-symbiont binary association matrix
  row_names <- cluster_sub$sample_id
  col_names <- cluster_sub$original_bin
  col_prefixes <- sub("_.*", "", col_names)
  assoc_mat <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                      dimnames = list(row_names, col_names))
  for (i in seq_along(row_names)) {
    for (j in seq_along(col_names)) {
      if (row_names[i] == col_prefixes[j]) {
        assoc_mat[i, j] <- 1
      }
    }
  }
  assoc_df <- as.data.frame(assoc_mat)
  assoc_df$sample_id <- rownames(assoc_df)
  assoc_df <- merge(assoc_df, input_meta_rtapa, by = "sample_id", all.x = TRUE)
  binary_assoc <- assoc_df %>%
    select(-sample_id) %>%
    group_by(population) %>%
    summarise(across(everything(), ~ as.integer(any(. > 0))))
  assoc_matrix <- as.data.frame(binary_assoc)
  colnames(assoc_matrix) <- sub("\\..*", "", colnames(assoc_matrix))
  rownames(assoc_matrix) <- assoc_matrix$population
  assoc_matrix$population <- NULL
  assoc_matrix <- as.matrix(assoc_matrix)
  
  # Filter the host tree to match taxa in the matrix
  tips_to_keep <- intersect(tree$tip.label, rownames(assoc_matrix))
  tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
  
  # Compute distance matrices
  gdist <- cophenetic(tree_filtered)
  ldist <- cophenetic(tree_strain)
  D <- prepare_paco_data(gdist, ldist, assoc_matrix)
  D <- add_pcoord(D)
  D <- PACo(D, nperm = 1000, seed = 24, symmetric = FALSE, method = "r0")
  
  # Extract and store PACo global statistics
  gof <- as.data.frame(D$gof)
  gof$cluster <- cluster_id
  global_results <- bind_rows(global_results, gof)
  
  # Run jackknife to evaluate individual interactions
  jackknife <- paco_links(D)
  dir.create(paste0("Results/Paco/", cluster_id), showWarnings = FALSE)
  write.csv(jackknife$jackknife,
            file = paste0("Results/Paco/", cluster_id, "/jackknife_scores.csv"),
            quote = FALSE)
  
  message("Cluster ", cluster_id, " processed successfully.")
}

# Save PACo global summary results
colnames(global_results) <- c("Paco_p", "Paco_ss", "Paco_n", "cluster")
write.csv(global_results,
          file = "Results/Paco/paco_summary_results.csv",
          row.names = FALSE)

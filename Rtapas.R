################## Rtapas ######################


########################
## Define Input Files ##
########################
library(ape)
library(dplyr)
library(tidyr)
library(Rtapas)
library(stringr)

df <- read.delim("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/data/Table_final_OTU.txt", header = TRUE, sep= "\t")
meta_rtapas <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ";")
cluster_tree <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Tree/Cdb.csv")
tree_strain_test <- read.tree(file = "data/Phylogeny_geography/Alltrees/1204_1_mugsy_tree.treefile")

# Chaîne Newick correcte
tree_human <- "(Baka,(Hadza,((California,Norman),(OrangAsli,(UrbanMalaysians,(Nepali,(Matses, Yanomami)))))));"

# Lire l’arbre
tree <- read.tree(text = tree_human)
tree <- compute.brlen(tree, method = 1)
# L’arbre est ladderisé pour que la racine soit en bas à gauche
tree <- ladderize(tree, right = FALSE)

# Affichage
plot(tree, main = "Population tree", direction = "rightwards", cex = 0.9)

# separate by cluster 
cluster_ratapas<-read.csv( "data/cluster_strains_final.csv ", header = T)
cluster_rtapas <- cluster_ratapas %>% 
  select(original_bin,secondary_cluster) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>% 
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  mutate(sample_id = original_bin)%>%
  mutate(sample_id = str_remove(sample_id, "_.*$"))

#Filter les datas
rtapas_filtré <- cluster_rtapas %>%
  filter(secondary_cluster == "1204_1")


# transform the table to be in long format 
# a tab-delimited file with samples as rows and metadata as columns

##################
## Filter Reads ##
##################
reads <- df[,sampleid]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(comm.df = reads)
dim(reads)

input.rtapas <- as.data.frame(t(reads))
#input.rtapas <- input.comm[sampleid,]
input.rtapas$sample_id <- rownames(input.rtapas)

rownames(meta_rtapas) <-meta_rtapas$sample_id
input_meta_rtapa <- meta_rtapas %>%
  select(sample_id,population) %>%
  mutate(population = case_when(
    population == "Urban malaysian" ~ "UrbanMalaysians",
    population == "Orang Asli" ~ "OrangAsli",
    TRUE ~ population
  ))
  
##merging dataset   
merged_rtapa <- merge(input.rtapas, input_meta_rtapa, by = "sample_id")


library(dplyr)

# creation empty matrice 

empty_df <- matrix(NA,
                   nrow = length(rtapas_filtré$sample_id),
                   ncol = length(rtapas_filtré$original_bin),
                   dimnames = list(rtapas_filtré$sample_id, rtapas_filtré$original_bin)) 
write.csv(empty_df, "Results/Rtapas/matrice.csv", na="NA")
###### Creation df to do the matrix of association 
# extraction of colonms names and row names 
row_names <- rownames(empty_df)
col_names <- colnames(empty_df)

# Extraction of the prefix before _. 
col_prefixes <- sub("_.*", "", col_names)

# Intiaitialise matrice to 0 
empty_df[,] <- 0

# fill with 1 if col_prefixes= row_names
for (i in seq_along(row_names)) {
  for (j in seq_along(col_names)) {
    if (row_names[i] == col_prefixes[j]) {
      empty_df[i, j] <- 1
    }
  }
}
#creation of the dataset with population 
test <- as.data.frame(empty_df)
test$sample_id <- rownames(test)

test <- merge(test, input_meta_rtapa[, c("sample_id", "population")],
             by= "sample_id", all.x = TRUE)


#group by pop and determinetest
binary_assoc <- test %>%
  select(-sample_id) %>%
  group_by(population) %>%
  summarise(across(everything(), ~ as.integer(any(. > 0)))) 


#create the matrix 
assoc_matrix <- as.data.frame(binary_assoc)
colnames(assoc_matrix)<- sub(pattern = "\\..*", replacement = "", x = colnames(assoc_matrix))
rownames(assoc_matrix) <- assoc_matrix$population
assoc_matrix$population <- NULL
assoc_matrix <- as.matrix(assoc_matrix)

write.table(assoc_matrix, file = "Results/Rtapas/matrice.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.csv(assoc_matrix, file = "Results/Rtapas/matrice.csv", row.names = TRUE)

###Rtapas 
##### see number of n 

# sum(assoc_matrix)* 0.1
# sum(assoc_matrix)* 0.2
# n <- 4

###Keep the same tips and names everywhere 
tips_to_keep <- intersect(tree$tip.label, rownames(assoc_matrix))

tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))


# Tips in host and symbiont trees
host_tips <- tree_filtered$tip.label
symbiont_tips <- colnames(assoc_matrix)

# Get number of unique one-to-one possible associations
# i.e., one host per symbiont (and vice versa)
# Compute n
max_n_possible <- min(length(tree_filtered$tip.label), length(tree_strain_test$tip.label))
n_prop <- round(sum(assoc_matrix) * 0.15)
n_max <- min(n_prop, max_n_possible)

message(paste("Calculated n =", n_max, "(capped at", max_n_possible, ")"))


PACo_LFi <- max_incong(
             HS = assoc_matrix, 
             treeH = tree_filtered, 
             treeS = tree_strain_test, 
             n= 4, 
             N= 1000, 
             method = "paco", 
             symmetric = FALSE,
             ei.correct = "sqrt.D", 
             percentile = 0.99, 
             diff.fq = TRUE,
             strat = "parallel" ,
             cl = 10)

 head(PACo_LFc)

 # Compute G*
G_star <- gini_RSV(PACo_LFi$LFr)
PACo_LFi$LFr <-sort( x= PACo_LFi$LFr)
barplot(PACo_LFi$LFr)
write.table(PACo_LFc, "Results/Rtapas/table_rtapas.txt")

##########PLOT tangle gram 

colgrad <- c("darkred", "gray90", "darkblue")

#create plot 

pdf(file = "Results/Rtapas/tangle_gram.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

tangle_gram(treeH = tree_filtered, treeS = tree_strain_test, HS = assoc_matrix, PACo_LFi, colscale = "diverging", 
            colgrad = c("darkred", "gray90", "darkblue"), nbreaks = 50, 
            node.tag = TRUE, cexpt = 1.2, link.lwd = 1, link.lty = 1, 
            fsize = 0.5, pts = FALSE, ftype ="i")


dev.off()

#####################

########################
## Define Input Files ##
########################
library(ape)
library(dplyr)
library(tidyr)
library(Rtapas)
library(stringr)
library(paco)

df <- read.delim("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/data/Table_final_OTU.txt", header = TRUE, sep= "\t")
meta_rtapas <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ";")
cluster_tree <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Tree/Cdb.csv")
tree_strain_test <- read.tree(file = "data/Phylogeny_geography/Alltrees/506_1_mugsy_tree.treefile")

# Chaîne Newick correcte
tree_human <- "(Baka,(Hadza,((California,Norman),(OrangAsli,(UrbanMalaysians,(Nepali,(Matses, Yanomami)))))));"

# Lire l’arbre
tree <- read.tree(text = tree_human)
tree <- compute.brlen(tree, method = 1)
# L’arbre est ladderisé pour que la racine soit en bas à gauche
tree <- ladderize(tree, right = FALSE)

# Affichage
plot(tree, main = "Population tree", direction = "rightwards", cex = 0.9)

# separate by cluster 
cluster_ratapas<-read.csv( "data/cluster_strains_final.csv ", header = T)
cluster_rtapas <- cluster_ratapas %>% 
  select(original_bin,secondary_cluster) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>% 
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  mutate(sample_id = original_bin)%>%
  mutate(sample_id = str_remove(sample_id, "_.*$"))

#Filter les datas
rtapas_filtré <- cluster_rtapas %>%
  filter(secondary_cluster == "506_1")


# transform the table to be in long format 
# a tab-delimited file with samples as rows and metadata as columns

##################
## Filter Reads ##
##################
reads <- df[,sampleid]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(comm.df = reads)
dim(reads)

input.rtapas <- as.data.frame(t(reads))
#input.rtapas <- input.comm[sampleid,]
input.rtapas$sample_id <- rownames(input.rtapas)

rownames(meta_rtapas) <-meta_rtapas$sample_id
input_meta_rtapa <- meta_rtapas %>%
  select(sample_id,population) %>%
  mutate(population = case_when(
    population == "Urban malaysian" ~ "UrbanMalaysians",
    population == "Orang Asli" ~ "OrangAsli",
    TRUE ~ population
  ))

##merging dataset   
merged_rtapa <- merge(input.rtapas, input_meta_rtapa, by = "sample_id")


library(dplyr)

# creation empty matrice 

empty_df <- matrix(NA,
                   nrow = length(rtapas_filtré$sample_id),
                   ncol = length(rtapas_filtré$original_bin),
                   dimnames = list(rtapas_filtré$sample_id, rtapas_filtré$original_bin)) 

###### Creation df to do the matrix of association 
# extraction of colonms names and row names 
row_names <- rownames(empty_df)
col_names <- colnames(empty_df)

# Extraction of the prefix before _. 
col_prefixes <- sub("_.*", "", col_names)

# Intiaitialise matrice to 0 
empty_df[,] <- 0

# fill with 1 if col_prefixes= row_names
for (i in seq_along(row_names)) {
  for (j in seq_along(col_names)) {
    if (row_names[i] == col_prefixes[j]) {
      empty_df[i, j] <- 1
    }
  }
}
#creation of the dataset with population 
test <- as.data.frame(empty_df)
test$sample_id <- rownames(test)

test <- merge(test, input_meta_rtapa[, c("sample_id", "population")],
              by= "sample_id", all.x = TRUE)


#group by pop and determinetest
binary_assoc <- test %>%
  select(-sample_id) %>%
  group_by(population) %>%
  summarise(across(everything(), ~ as.integer(any(. > 0)))) 


#create the matrix 
assoc_matrix <- as.data.frame(binary_assoc)
colnames(assoc_matrix)<- sub(pattern = "\\..*", replacement = "", x = colnames(assoc_matrix))
rownames(assoc_matrix) <- assoc_matrix$population
assoc_matrix$population <- NULL
assoc_matrix <- as.matrix(assoc_matrix)



###paco

###Keep the same tips and names everywhere 
tips_to_keep <- intersect(tree$tip.label, rownames(assoc_matrix))

tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))


# Tips in host and symbiont trees
host_tips <- tree_filtered$tip.label
symbiont_tips <- colnames(assoc_matrix)

gdist <- cophenetic(tree_filtered)
ldist <- cophenetic(tree_strain_test)
D <- prepare_paco_data(gdist, ldist, assoc_matrix)
D <- add_pcoord(D)
D <- PACo(D, nperm=1000, seed=42, method="backtrack", symmetric = FALSE)
print(D$gof)

residus <- residuals_paco(D$proc,type = "interaction")
hist(residus, main = "Distribution des résidus Procrustes", xlab = "Résidu", col = "gray70", border = "white")



jacknife <-paco_links(D, .parallel = FALSE, proc.warnings = TRUE)
jackknife_scores <- jacknife$jackknife
hist(jackknife_scores, main = "Jackknife scores (interaction contribution)", col = "lightblue", xlab = "Score")

    

result_paco <- as.data.frame(D$gof)

######

library(ape)
library(dplyr)
library(tidyr)
library(Rtapas)
library(stringr)

# Load abundance and metadata
df <- read.delim("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/data/Table_final_OTU.txt", header = TRUE, sep= "\t")
meta_rtapas <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ";")
cluster_ratapas <- read.csv("data/cluster_strains_final.csv ", header = TRUE)

# Format cluster table
cluster_rtapas <- cluster_ratapas %>%
  select(original_bin, secondary_cluster) %>%
  mutate(across(everything(), ~ sub("^.*__", "", .))) %>%
  mutate(across(secondary_cluster, ~ gsub("\\.", "", .))) %>%
  mutate(sample_id = str_remove(original_bin, "_.*$"))

# Human tree
tree_human <- "(Baka,(Hadza,((California,Norman),(OrangAsli,(UrbanMalaysians,(Nepali,(Matses, Yanomami)))))));"
tree <- read.tree(text = tree_human)
tree <- ladderize(compute.brlen(tree, method = 1), right = FALSE)

# Reads CPM
reads <- df[, sampleid]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(comm.df = reads)

# Input matrix for Rtapas
input.rtapas <- as.data.frame(t(reads))
input.rtapas$sample_id <- rownames(input.rtapas)

# Fix metadata
rownames(meta_rtapas) <- meta_rtapas$sample_id
input_meta_rtapa <- meta_rtapas %>%
  select(sample_id, population) %>%
  mutate(population = case_when(
    population == "Urban malaysian" ~ "UrbanMalaysians",
    population == "Orang Asli" ~ "OrangAsli",
    TRUE ~ population
  ))

# Merge abundance and metadata
merged_rtapa <- merge(input.rtapas, input_meta_rtapa, by = "sample_id")

# List of all tree files
tree_files <- list.files("data/Phylogeny_geography/Alltrees", pattern = "_mugsy_tree.treefile$", full.names = TRUE)

for (tree_path in tree_files) {
  try({
    cluster_id <- gsub("_mugsy_tree.treefile", "", basename(tree_path))
    message(" Processing cluster: ", cluster_id)
    
    # Filter cluster
    rtapas_filtré <- cluster_rtapas %>%
      filter(secondary_cluster == cluster_id)
    
    if (nrow(rtapas_filtré) == 0) next
    
    # Build binary association matrix
    empty_df <- matrix(0,
                       nrow = length(rtapas_filtré$sample_id),
                       ncol = length(rtapas_filtré$original_bin),
                       dimnames = list(rtapas_filtré$sample_id, rtapas_filtré$original_bin))
    
    row_names <- rownames(empty_df)
    col_names <- colnames(empty_df)
    col_prefixes <- sub("_.*", "", col_names)
    
    for (i in seq_along(row_names)) {
      for (j in seq_along(col_names)) {
        if (row_names[i] == col_prefixes[j]) {
          empty_df[i, j] <- 1
        }
      }
    }
    
    # Add population info
    test <- as.data.frame(empty_df)
    test$sample_id <- rownames(test)
    
    test <- merge(test, input_meta_rtapa[, c("sample_id", "population")],
                  by = "sample_id", all.x = TRUE)
    
    binary_assoc <- test %>%
      select(-sample_id) %>%
      group_by(population) %>%
      summarise(across(everything(), ~ as.integer(any(. > 0))))
    
    assoc_matrix <- as.data.frame(binary_assoc)
    colnames(assoc_matrix) <- sub(pattern = "\\..*", replacement = "", x = colnames(assoc_matrix))
    rownames(assoc_matrix) <- assoc_matrix$population
    assoc_matrix$population <- NULL
    assoc_matrix <- as.matrix(assoc_matrix)
    
    # Tree processing
    tips_to_keep <- intersect(tree$tip.label, rownames(assoc_matrix))
    tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
    tree_strain <- read.tree(tree_path)
    
    # Compute dynamic n
    n_prop <- round(sum(assoc_matrix) * 0.15)
    n_max <- min(length(tree_filtered$tip.label), length(tree_strain_test$tip.label))
    n <- min(n_prop, n_max)
    
    message(" Calculated n = ", n)
    
    
    # Run Rtapas
    PACo_LFi <- max_incong(
      HS = assoc_matrix, 
      treeH = tree_filtered, 
      treeS = tree_strain, 
      n= 3, 
      N= 8, 
      method = "paco", 
      symmetric = TRUE,
      ei.correct = "sqrt.D", 
      percentile = 0.99, 
      diff.fq = TRUE, 
      strat = "parallel",
      cl = 10)
    
    PACo_LFc <- max_cong(
      HS = assoc_matrix,
      treeH = tree_filtered,
      treeS = tree_strain,
      n = 3,
      N = 8,
      method = "paco",
      symmetric = TRUE,
      ei.correct = "sqrt.D",
      percentile = 0.01,
      res.fq = TRUE,
      strat = "parallel",
      cl = 10
    )
    
    
    # Create output folder
    output_dir <- file.path("Results/Rtapas", cluster_id)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    write.table(PACo_LFi, file = file.path(output_dir, "table_rtapas.txt"))
    
    # Plot tanglegram
    pdf(file = file.path(output_dir, "tangle_gram.pdf"),
        width = 10, height = 10)
    
    tangle_gram(treeH = tree_filtered,
                treeS = tree_strain,
                HS = assoc_matrix,
                fqtab = PACo_LFi,
                colscale = "sequential",
                colgrad = c("darkred", "gray90", "darkblue"),
                nbreaks = 50,
                node.tag = TRUE,
                cexpt = 1.2,
                link.lwd = 1,
                link.lty = 1,
                fsize = 0.5,
                pts = FALSE,
                ftype = "i")
    
    dev.off()
    message(" tapas completed for cluster ", cluster_id)
  }, silent = TRUE)
}






# Test decreasing values of N until no warning
# run_max_cong_variable_N <- function(assoc_matrix, treeH, treeS, n, N_values = c(3000,2000,1000, 500, 200,150, 100,90,80,70,60, 50,40,30,20,15,10)) {
#   for (N in N_values) {
#     warning_triggered <- FALSE
#     
#     PACo_LFc <- withCallingHandlers(
#       {
#         res <- max_cong(
#           HS = assoc_matrix,
#           treeH = treeH,
#           treeS = treeS,
#           n = n,
#           N = N,
#           method = "paco",
#           symmetric = TRUE,
#           ei.correct = "sqrt.D",
#           percentile = 0.01,
#           res.fq = TRUE,
#           cl = 4
#         )
#         if (!warning_triggered) {
#           message(paste("Rtapas successful with N =", N, "and n =", n))
#           return(list(result = res, N = N, warning = FALSE))
#         }
#       },
#       warning = function(w) {
#         if (grepl("No. of trimmed H-S assoc. matrices", conditionMessage(w))) {
#           warning_triggered <<- TRUE
#           message(paste("️ Warning with N =", N, ":", conditionMessage(w)))
#           invokeRestart("muffleWarning")  # skip printing warning
#         }
#       }
#     )
#   }
#   stop(paste("Rtapas failed for all N values with n =", n))
# }


# count_one_to_one_options <- function(assoc_matrix, n) {
#   hosts <- rownames(assoc_matrix)
#   symbs <- colnames(assoc_matrix)
#   
#   valid_pairs <- which(assoc_matrix == 1, arr.ind = TRUE)
#   unique_hosts <- length(unique(rownames(assoc_matrix)[valid_pairs[, 1]]))
#   unique_symbs <- length(unique(colnames(assoc_matrix)[valid_pairs[, 2]]))
#   
#   possible <- min(unique_hosts, unique_symbs)
#   return(possible >= n)
# }
# 
# 
# run_max_cong_variable_N(assoc_matrix = assoc_matrix, treeH = tree_filtered, treeS = tree_strain_test, n=2)
# library(pheatmap)
# pheatmap(assoc_matrix, cluster_rows = FALSE, cluster_cols = FALSE, main = "Binary association matrix")

#PACo_LFc <- max_cong(HS = assoc_matrix, treeH = tree_filtered, treeS= tree_strain_test, n = 3, N = 10, method = "paco", symmetric = TRUE, 
#  ei.correct = "sqrt.D", percentile = 0.01, res.fq = TRUE, 
#     cl = 10)

# PACo_LFc <- max_cong(
#   HS = assoc_matrix,
#   treeH = tree_filtered,
#   treeS = tree_strain_test,
#   n = 5,
#   N = 500,
#   method = "paco",
#   symmetric = TRUE,
#   ei.correct = "sqrt.D",
#   percentile = 0.01,
#   res.fq = TRUE,
#   cl=10
# )

######Maaslin2
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")
rm(list=ls())
library(Maaslin2)
####################
# custom functions #
####################
convert_to_CPM <- function(comm.df){
  comm.df <- comm.df*10^6
  comm.df}
prevalence_filter <- function(comm.df, n) {
  x <- comm.df # matrix of abundance with sampleid as columns
  y <- x
  y[y>0] <- 1
  y <- y[rowSums(y)>n,]
  x <- x[rownames(x) %in% rownames(y),]
  x
}
headx <- function(data.frame,n=NULL){
  if(is.null(n)){
    data.frame[1:5,1:5]
  }else{
    data.frame[1:n,1:n]
  }
}

#load data 
df <- read.delim("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/data/Table_final_OTU.txt", header = TRUE, sep= "\t")
def_meta <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ",")
ls.var <- c("population", "Age_bin","Age","Gender","Country","Lat.Long","Village","Lifestyle","BMI","State")

###################
## Select Sample ##
###################
samples_exclude <- c("ERR7738433", "ERR7738234", "ERR7738665") # ERR7738433_1.fastq.gz and ERR7738234_2.fastq.gz were not downloaded successfully, hence all downstream steps with this sample failed. ERR7738665 paired reads have different names.
sampleid <- def_meta$sample_id [! def_meta$sample_id %in% samples_exclude]
metadata <- def_meta[def_meta$sampleid %in% colnames(df),]
metadata <- def_meta
#bio <- bio[bio$sampleid %in% colnames(df),c("sampleid","Age","Gender")]
#ffq <- ffq[ffq$sampleid %in% colnames(df),c("sampleid",ls.var)]

#metadata <- merge.data.frame(metadata, bio, by="sampleid")
#metadata <- merge.data.frame(metadata, ffq, by="sampleid")

#Ob <- bio$sampleid[bio$obesity=="obese"&!is.na(bio$obesity)]
#incomplete <- as.character(ffq$sampleid)[apply(ffq,1, function(x) any(is.na(x)))]

#sampleid <- as.character(metadata$sampleid)
#sampleid <- sampleid[!sampleid %in% incomplete]
#sampleid <- sampleid[!sampleid %in% Ob]

##################
## Filter Reads ##
##################
reads <- df[,sampleid]
rownames(reads) <- df$Genome
reads <- convert_to_CPM(comm.df = reads)
reads <- prevalence_filter(comm.df = reads, n = 5) # abundance >0, prevalence >5
dim(reads)

########################
## Define Input Files ##
########################
# a tab-delimited file with samples as rows and metadata as columns

input.comm <- as.data.frame(t(reads))
input.comm <- input.comm[sampleid,]
input.metadata <- metadata[metadata$sample_id %in% sampleid, ]
rownames(input.metadata) <-input.metadata$sample_id
input.metadata <- input.metadata[sampleid,]

###############
## Fit Model ##
###############
library(Maaslin2)
print(paste("Using", paste(ls.var, collapse = ", "), "as fixed effect (ref. PERMANOVA Model Selection)."))

# set factor and level
# the first item is automatically set as reference
pop.ord <- as.character(unique(input.metadata$population))
n_basap <- which(pop.ord=="Hadza")
pop.ord <- c(pop.ord[n_basap],pop.ord[-n_basap])
input.metadata$population <- factor(input.metadata$population, pop.ord)
input.metadata$Age_bin <- factor(input.metadata$Age_bin, levels = c("18 and above", "14 and below"))
input.metadata$Country <- factor(input.metadata$Country, levels =c("Tanzania","Nepal","USA","Brazil","Cameroon","Peru","Malaysia"))
input.metadata$Author <- factor(input.metadata$Author, levels =c("Carter","Conteville","Rampelli","Obregon-Tito","Tee"))
input.metadata$Lifestyle <- factor(input.metadata$Lifestyle, levels =c("H-G","Industrial"))
input.metadata$Geographic_region <- factor(input.metadata$Lifestyle, levels =c("Africa","Asia","North_Amercia","South_America"))
print("Hadza was used as the reference population.")
print("18 and above was used as the reference age.")
write.table(input.comm, "data//maaslin_input_comm.tsv",quote = F,sep = "\t")
write.table(input.metadata,"data//maaslin_input_metadata.tsv",quote = F,sep = "\t", row.names = FALSE)


#########################
# fit univariate models #
#########################

# population as fixed effect
i="Geographic_region"
output_dir0 <- "Results/Maaslin/"
output_dir <- paste0(output_dir0,i,"_author_pop_ranef")

fit <- Maaslin2(input_data = input.comm,
                input_metadata = input.metadata,
                output=output_dir,
                fixed_effects = c("Geographic_region"),
                normalization = "CLR",
                transform = "None",
                plot_scatter = F,
                save_models=T,
                random_effects = c("Author", "population"),
                cores = 3)

# Lifestyle as fixed effect
i="N_reads_clean"
output_dir0 <- "Results/Maaslin/"
output_dir <- paste0(output_dir0,i,"_pop_ranef")

fit <- Maaslin2(input_data = input.comm,
                input_metadata = input.metadata,
                output=output_dir,
                fixed_effects = c("N_reads_clean"),
                normalization = "CLR",
                transform = "None",
                plot_scatter = F,
                save_models=T,
                random_effects = c("population"),
                cores = 3)

ls.var <- c( "Age_bin","Gender","Country","Lifestyle")

for(i in ls.var){
  output_dir_wpop <- paste0(output_dir0,i,"_pop_ranef")
  
  print(paste("Started MAASLIN2 for", i))
  
  
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_wpop,
                  fixed_effects = c(i),
                  random_effects = c("population"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = T,
                  save_models=T,
                  cores = 3)
  
  print(paste(i, "finished running."))
}


# prelim mag count
signif_count <- data.frame()
signif_list <- list()
for(i in c(ls.var)){
  output_dir0 <- "Results/Maaslin/"
  for(m in c("_pop_ranef")){
    input.file <- paste0(output_dir0,i,m,"/all_results.tsv")
    input.file <- read.delim(input.file)
    input.file <- input.file[input.file$qval<0.05, ]
    signif.features <- unique(input.file$feature)
    
    
    if(length(signif.features)==0){
      signif.features <- c("No Associations")
      out <- data.frame(model=paste0(i,m), var=i, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=0)
      signif_count <- rbind(signif_count, out)
    } else {
      out <- data.frame(model=paste0(i,m), var=i, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=length(signif.features))
      signif_count <- rbind(signif_count, out)
    }
    
    len <- length(signif_list)
    
    if(len==0){
      signif_list[[1]] <- signif.features
      names(signif_list) <- paste0(i,m)
    } else {
      nm <- c(names(signif_list),paste0(i,m))
      signif_list[[len+1]] <- signif.features
      names(signif_list) <- nm
    }
  }
}


##############################
# select multivariate models #
##############################
multivar_list <- unique(signif_count[signif_count$signif.feature!=0,"var"])
#exclude_var <- c("Hunting", "Watching_TV","population")
#multivar_list <- multivar_list[!multivar_list%in% exclude_var]

AICcPermanova::make_models(multivar_list,k = 3) -> my.models

tmp <- unlist(lapply(strsplit(my.models$form ,split = "~"), function(x) x[2]))
model.list <- lapply(strsplit(tmp, split="\\+"), function(x) gsub(pattern = " ", replacement = "", x))

model.list[which(lapply(model.list, function(x) length(x)) >1)] -> model.list

output_dir0 <- "Results/Maaslin/"
for(i in 1:length(model.list)){
  ivar <- model.list[[i]]
  pref <- paste0(ivar,collapse = "+")
  #output_dir <- paste0(output_dir0,pref,"_noranef")
  output_dir_wpop <- paste0(output_dir0,pref,"_pop_author_ranef")
  
  print(paste("Started MAASLIN2 for", pref))
  
  
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_wpop,
                  fixed_effects = c(ivar),
                  random_effects = c("population", "Author"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = T,
                  save_models=T,
                  cores = 3)
  
  print(paste(pref, "finished running."))
}


for(i in 1:length(model.list)){
  ivar <- model.list[[i]]
  pref <- paste0(ivar,collapse = "+")
  
  output_dir0 <- "Results/Maaslin/"
  
  for(m in c("_pop_ranef")){
    input.file <- paste0(output_dir0,pref,m,"/all_results.tsv")
    input.file <- read.delim(input.file)
    input.file <- input.file[input.file$qval<0.05, ]
    signif.features <- unique(input.file$feature)
    
    if(length(signif.features)==0){
      signif.features <- c("No Associations")
      out <- data.frame(model=paste0(pref,m), var=pref, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=0)
      signif_count <- rbind(signif_count, out)
    } else {
      out <- data.frame(model=paste0(pref,m), var=pref, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=length(signif.features))
      signif_count <- rbind(signif_count, out)
    }
    
    len <- length(signif_list)
    nm <- c(names(signif_list),paste0(pref,m))
    signif_list[[len+1]] <- signif.features
    names(signif_list) <- nm
  }
}



str(signif_list)
signif_count[signif_count$signif.feature>0,]

saveRDS(signif_list, file = "Results/Maaslin/significant.MAGs.list.rds")
MAG.list <- signif_list
saveRDS(signif_list, file = "Results/Maaslin/significant.MAGs.list.rds")



#######################
## MAASLIN author ranef ##
#######################
library(lme4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load model 
models <- readRDS("Results/Maaslin/population_author_ranef/fits/models.rds")
pop.varex <- data.frame(feature=c(), pop_var=c())
for (i in names(models)){
  test <- models[[i]]
  ranef.var <- data.frame(VarCorr(test))[,c("grp","var1","vcov","sdcor")]
  varex <- ranef.var$vcovvqr
  names(varex) <- ranef.var$grp
  varex <- round(varex*100/sum(varex),1)
  
  tmp <- data.frame(feature=i,pop_var=varex[1],row.names = NULL)
  
  pop.varex <- rbind(pop.varex,tmp)
}



# most abundant in
#df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund <- input.comm
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

#input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata_2<- input.metadata
input.metadata_2$sampleid <- rownames(input.metadata_2)
input.metadata_2 <- input.metadata_2[,c("sampleid","population")]
full.taxonomy <- Abund_strong_2 %>%
  select("Genome", "phylum","class","order","family","genus","species") %>%
  rename(feature = Genome)
df.abund <- merge.data.frame(input.metadata_2,df.abund, by="sampleid")

df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")



df.abund %>% group_by(phylum, feature, population) %>% summarise(popMean=mean(CPM)) -> pop_average

mostabund.pop <- data.frame(feature=c(), topPop=c(), meanCPM=c())
for(i in unique(pop_average$feature)){
  tmp <- pop_average[pop_average$feature==i,]
  max <- max(tmp$popMean)
  n <- grep(max, tmp$popMean)
  
  j <- tmp$population[n]
  
  tmp <- data.frame(feature=i, topPop=j, meanCPM=max)
  mostabund.pop <- rbind(mostabund.pop, tmp)
}
mostabund.pop

popranef.df <- merge.data.frame(pop.varex, mostabund.pop, by="feature")
write.table(popranef.df, 
            "Results/Maaslin/population_author_ranef/all_results_popvar_author_ranef.tsv", 
            quote = F,sep = "\t", row.names = F)

popvar.cutoff <- 35 # in percentage
strong.popranef <- unique(popranef.df[popranef.df$pop_var>=popvar.cutoff,"feature"])
length(strong.popranef)


##################
## diet volcano ##
##################
df <- read.delim("Results/Maaslin/population_author_ranef/all_results.tsv")
#df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_noranef/all_results_wTax.tsv")
length(unique(df$feature)) -> ct
df <- subset(df, qval < 0.05)
length(unique(df$feature)) -> sig.ct

sig.ct *100 / ct
df$Genome <- df$feature
df <- merge(df, Phylum, by= "feature" )

ggdata <- df
ggdata$association <- ggdata$value
ggdata$association[ggdata$qval>=0.05] <- "none"
ggdata$association <- factor(ggdata$association, c("Orang Asli","California","Nepali","Baka","Matses","Yanomami","Norman","Urban malaysian","none"))

write.table(ggdata ,"Results/Maaslin/population_author_ranef/allresults_info.txt",
            quote = F, row.names = F, sep = "\t")
#ggdata <- df[,c("feature","value","coef","N.not.0","pval","qval","species","Genome")]
ggdata <- ggdata[order(abs(ggdata$coef),decreasing = T),]
# ggdata$label <- seq(1,nrow(ggdata),by=1)
ggdata$label <- ggdata$species
ggdata$label[ggdata$qval>=0.05] <- ""
rownames(ggdata) <- NULL


diet.label <- ggdata[ggdata$label!="",]
diet.label <- merge.data.frame(diet.label, full.taxonomy[,c("feature","phylum","family")], by="feature")
rownames(diet.label) <- NULL
write.table(diet.label,"Results/Maaslin/population_author_ranef/all_results_pop_label.tsv",
            quote = F, row.names = F, sep = "\t")

volc <- ggplot(ggdata, aes(y=-log10(qval), x=coef, fill=association)) + theme_bw() +
  geom_text(aes(label=label), col="grey70", size=3,nudge_y =0.1) +
  annotate("text", x = 0, y = -log10(0.045), label = "qval=0.05", vjust=0, col="red")+
  geom_hline(yintercept = -log10(0.05), col="red", lty=2)+
  geom_point(size=5, pch=21, col="grey30") +
  scale_fill_manual(values = c("#2E8B57","#FF69B4","#63B8FF","grey70","#20B2AA","#8B008B","#FF7F50","#DC143C","#FFD700")) +
  scale_x_continuous(limits = c(-3,3), breaks = sort(c(seq(-4,4, by=1),0)))+
  scale_y_continuous(limits=c(0, 10)) +
  labs(title = "Diet Correlations (Multi-variate)",
       caption = "q-val is p-val adjusted with Benjamin Hochberg method, q-val<0.05 considered as significant") +
  ylab("Significance -log10(qval)") + xlab("Coefficient")

volc

outpath <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/fig_out/"
for(d in c("png","pdf")){
  ggsave(plot = volc, 
         filename = paste0(outpath,"diet_volcano.",d), scale = 1.5,
         device = d,width = 10,height = 6,units = "in",dpi = "print")}

# Compare counts
#significant.count <- readRDS("Results/Maaslin/population_author_ranef/significant_results.tsv")
significant.count <- read.delim("Results/Maaslin/population_author_ranef/significant_results.tsv")
library(ggplot2)
library(ggpubr)
library(ggvenn)

############################################################

library(ggplot2)
library(ggpubr)



# Nettoyer les noms de colonnes si besoin
colnames(significant.count) <- gsub('"', '', colnames(significant.count))

# Sélectionner seulement les MAGs significatifs (qval < 0.05)
significant.features <- significant.count[significant.count$qval < 0.05, ]

# Compter le nombre de MAGs par population (value)
feature_count <- as.data.frame(table(significant.features$value))
colnames(feature_count) <- c("population", "signif.feature")

# Plot
countbar <- ggplot(feature_count) +
  theme_pubclean(base_size = 16) +
  geom_bar(aes(x = population, y = signif.feature), stat = "identity", fill = "#7A67EE", color = "black", width = 0.8) +
  geom_text(aes(x = population, y = signif.feature + 2, label = signif.feature)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Count of Significant MAGs per Population", subtitle = "q-value < 0.05", y = "Number of Significant MAGs", x = "")

# Afficher le plot
print(countbar)

# Sauvegarder si besoin
outpath <- "Results/Maaslin/population_author_ranef/"
for(d in c("png","pdf")){
  ggsave(plot = countbar, 
         filename = paste0(outpath,"signif_count_bar.",d), scale = 1.5,
         device = d,width = 10,height = 6,units = "in",dpi = "print")}










###############
## diet list significant associations  ##
###############
df <- read.delim("Results/Maaslin/population_author_ranef/allresults_info.txt", header=T,sep= "\t")
ggdata<- df
#ggdata <- df[,c("feature","value","coef","pval","N.not.0","qval","phylum","class","family","species")]
ggdata <- ggdata[ggdata$qval<0.00000000005,]
diet.MAGs <- unique(ggdata$feature)
length(diet.MAGs)

# order data
ggdata <- ggdata[order(ggdata$coef,decreasing = F),]
ggdata$value <- factor(ggdata$value, c("Orang Asli","California","Nepali","Baka","Matses","Yanomami","Norman","Urban malaysian","none"))

mag.levels <- as.character(unique(ggdata$feature))
ggdata$feature <- factor(ggdata$feature, mag.levels)

#phylum.tree.order <- rev(readRDS("output_files/revised_MAG/phylum.tree.order.rds"))
phylum.tree.order <-unique(full.taxonomy$phylum)
#phylum.tree.order <- c(phylum.tree.order,"Methanobacteriota")
ggdata$phylum <- factor(ggdata$phylum, phylum.tree.order)

diet.list.df <- ggdata 

########
# plot #
########

#ggplot(diet.list.df, aes(x=coef,y=feature))+ theme_minimal() +
# theme(plot.margin = margin(0,0,0,0),
#      panel.spacing = unit(0, "cm"),
#     legend.position = "right",
#    panel.background = element_rect(color = "black",fill=NA),
#   strip.text.y.left = element_text(angle=0),
#  strip.background = element_rect(fill="white",linewidth = 0),
#       text=element_text(size=12), axis.text=element_text(size=12),
#      axis.text.y = element_text(angle=0, hjust = 1, size=8),
#     # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#    axis.title.y= element_blank(),
#   panel.grid.major.x = element_line(colour="grey",linetype = 2)) + 
#  scale_color_manual(values = c("#00A600","#FF69B4","#0000FF","grey70")) +
# scale_x_continuous(breaks = seq(-3,4, by=0.5),limits = c(-3,4), expand = c(0.05,0))+
#xlab("Coefficients") +
#  labs(title="MAGs with Significant Diet Associations (qval<0.05)") +
# facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#geom_point(aes(col=value),size=3) +
#  geom_text(aes(label=species), angle=0, nudge_x = 0.1, 
#           size=3, col="grey40", hjust=0, vjust=0.2)

#
# abundance heatmap (total)
library(reshape2)
library(dplyr)
#df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund <- input.comm
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

#input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata_2$sample_id <- rownames(input.metadata)
input.metadata_2 <- input.metadata[,c("sample_id","population")]
input.metadata_2 <- input.metadata_2 %>%
  rename(sampleid = sample_id)
df.abund <- merge.data.frame(input.metadata_2,df.abund, by="sampleid")
df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

n <- length(unique(df.abund$sampleid))

####Calculer l’abondance moyenne normalisée pour chaque MAG (feature).###
#####C’est une façon de produire une abondance relative moyenne, cohérente avec d’autres analyses####
df.abund %>% group_by(feature, phylum) %>% summarise(abundance=sum(CPM)/(n*10^6)) -> total_abund

heat.df.total <- total_abund[total_abund$feature %in% mag.levels,]
heat.df.total$feature <- factor(heat.df.total$feature, mag.levels)
heat.df.total$phylum <- factor(heat.df.total$phylum, phylum.tree.order)


#ggplot(heat.df.total)+ theme_minimal() + 
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_tile(aes(x="A", y=feature,fill=log(abundance)))+
#  scale_fill_gradient(low = "grey90", high = "black",)+
#  theme(plot.margin = margin(0,0,0,0),
#        panel.spacing = unit(0, "cm"),
#       legend.position = "right",
#      panel.grid = element_blank(),
#     panel.border = element_rect(colour = "black", fill=NA),
#    strip.text.y.left = element_text(angle=0),
#   strip.background = element_rect(fill="white",linewidth = 0),
#  text=element_text(size=12), axis.text=element_text(size=12),
# axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "white"),
#       # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#      axis.text.y = element_text(size=8),
#     axis.title.y= element_blank()) +
#  labs(title="Total Abundance") +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0))+
#  xlab("")

# abundance heatmap (bypop)
# normalised average
# the normalisation is the percentage of the population average over the sum of the average value in each population by features
#df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
#df.abund$sampleid <- rownames(df.abund)
#df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

#input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
#input.metadata$sampleid <- rownames(input.metadata)
#input.metadata <- input.metadata[,c("sampleid","population")]

#df.abund <- merge.data.frame(input.metadata_2,df.abund, by="sampleid")
#df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")


# n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(population,feature, phylum) %>% 
  summarise(n=length(unique(sampleid)),
            popCPM=sum(CPM),
            popMedian=median(CPM),
            popMean=mean(CPM),
            popSD=sd(CPM)) -> pop_abund

####total_abund contient l’abondance moyenne normalisée de chaque MAG dans chaque phylum.####

df.abund %>% group_by(feature) %>% summarise(totalCPM=sum(CPM)) -> total_abund

heat.df.pop <- merge.data.frame(pop_abund, total_abund, by="feature")
#####Tu additionnes les moyennes par population pour chaque MAG.####
####Cela te donne un dénominateur pour normaliser les moyennes (voir ci-dessous).####
pop_abund %>% group_by(feature) %>% summarise(sum.popMean=sum(popMean)) -> sumpopMean
heat.df.pop <- merge.data.frame(heat.df.pop, sumpopMean, by="feature")

heat.df.pop <- heat.df.pop[order(heat.df.pop$feature),]
####normalised.abund : proportion (%) de l’abondance totale du MAG qui vient d’une population donnée.###
###Qui domine en abondance brute (popCPM)###
###Quelle proportion de l’abondance totale (tous échantillons confondus) d’un MAG provient d’une population donnée ?

heat.df.pop$normalised.abund <- 100*heat.df.pop$popCPM/heat.df.pop$totalCPM
####normalised.average : contribution (%) de la moyenne d’une population par rapport à la somme des moyennes toutes populations confondues.###
###Qui a la moyenne la plus élevée (popMean), en valeur relative.###
####ne population a peu d’échantillons, mais le MAG est très abondant en moyenne chez elle → ça ne se voit pas avec popCPM, mais ça saute aux yeux avec normalised.average.
###Tu veux comparer l’intensité moyenne d’un MAG dans chaque population, pas juste "combien de lectures".
####
heat.df.pop$normalised.average <- 100*heat.df.pop$popMean/heat.df.pop$sum.popMean

# order data
heat.df.pop <- heat.df.pop[heat.df.pop$feature %in% mag.levels,]
heat.df.pop$feature <- factor(heat.df.pop$feature,mag.levels)
heat.df.pop$phylum <- factor(heat.df.pop$phylum, phylum.tree.order)
heat.df.pop$population <- factor(heat.df.pop$population,
                                 c("Hadza","Orang Asli","California","Nepali","Baka",
                                   "Matses","Yanomami","Norman","Urban malaysian"))

# ggplot(heat.df.pop, aes(x=population, y=feature, fill= normalised.average)) + 
# facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#theme_minimal() + 
#theme(plot.margin = margin(0,0,0,0),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.title.y= element_blank()) +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0))+
#  geom_tile()+ scale_fill_gradient(low="white",high="black") +
#  labs(title="Normalised Average") +
#  xlab("")

########
# plot #
########
library(ggplot2)
diet.cols <- c("#2E8B57","#FF69B4","#63B8FF","grey70","#20B2AA","#8B008B","#FF7F50","#DC143C","#FFD700")
names(diet.cols) <- c("Hadza","Orang Asli","California","Nepali","Baka",
                      "Matses","Yanomami","Norman","Urban malaysian")

a <- ggplot(diet.list.df, aes(x=coef,y=feature))+ 
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
  scale_color_manual(values = diet.cols) +
  scale_x_continuous(breaks = seq(-10,10, by=1),limits = c(-10,10), expand = c(0.06,0.06))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Coefficients") +
  labs(title="MAGs with Significant Diet Associations (qval<0.05)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_point(aes(col=value),size=3) +
  geom_text(aes(label=species), angle=0, nudge_x = 0.1, 
            size=3, col="grey40", hjust=0, vjust=0.2)

b <- ggplot(heat.df.total)+ 
  theme_minimal() +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_tile(aes(x="Total Abundance", y=feature,fill=log(abundance)))+
  scale_fill_gradient(low = "grey90", high = "black",)+
  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        #strip.text.y.left = element_text(angle=0),
        #strip.background = element_rect(fill="white",linewidth = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=10, colour = "black",angle=60, hjust = 1),
        axis.ticks.x= element_line(colour = NA),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),
        axis.title.y= element_blank()) +
  labs(title="") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("")

c <- ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  theme_minimal() + 
  theme(plot.margin = margin(0,0,2,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid = element_blank(),
        #strip.text.y.left = element_text(angle=0),
        #strip.background = element_rect(fill="white",linewidth = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=10, colour = "black",angle=60, hjust = 1),
        axis.ticks.x= element_line(colour = NA),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),
        axis.title.y= element_blank()) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  geom_tile()+ scale_fill_gradient(low="grey95",high="black") +
  labs(title="") +
  xlab("")

library(gridExtra)
lay <- rbind(c(rep(1,33),
               rep(2,1),
               rep(3,5)))

# draw plot
diet.list <- grid.arrange(a,b,c,layout_matrix=lay)
ggsave(filename = "Figures/Mags_with_pop_associations_qval_0.00000000005.pdf",diet.list,device = "pdf",width = 8.5, height = 9.5,units = "in")
saveRDS(grid.arrange(a,b,c,layout_matrix=lay),
        file = "Figures/Mags_with_pop_associations.rds")







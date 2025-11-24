install.packages("devtools")
library("devtools")
#install_github("biobakery/maaslin3")
BiocManager::install("biobakery/maaslin3")
library(maaslin3)

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
def_meta <- read.csv("/Users/lucie/Documents/Cambridge/Analysis/Analysis_thesis/Sample_pop_meta/metadata_pop_all_final.csv", header = TRUE, sep= ";")
#ls.var <- c("population", "Age_bin","Age","Gender","Country","Lat.Long","Village","Lifestyle","BMI","State")

###################
## Select Sample ##
###################
samples_exclude <- c("ERR7738433", "ERR7738234", "ERR7738665") # ERR7738433_1.fastq.gz and ERR7738234_2.fastq.gz were not downloaded successfully, hence all downstream steps with this sample failed. ERR7738665 paired reads have different names.
sampleid <- def_meta$sample_id [! def_meta$sample_id %in% samples_exclude]
metadata <- def_meta[def_meta$sampleid %in% colnames(df),]
metadata <- def_meta

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
input.metadata <- metadata
###############
## Fit Model ##
###############
library(maaslin3)
print(paste("Using", paste(ls.var, collapse = ", "), "as fixed effect (ref. PERMANOVA Model Selection)."))
### option
input.metadata <- input.metadata %>%
  filter(!population %in% c("Urban malaysian", "California", "Norman"))
input.metadata <- input.metadata %>%
  filter(population %in% c("Urban malaysian", "California", "Norman"))
# set factor and level
# the first item is automatically set as reference
pop.ord <- as.character(unique(input.metadata$population))
#n_basap <- which(pop.ord=="California")
n_basap <- which(pop.ord=="Hadza")
pop.ord <- c(pop.ord[n_basap],pop.ord[-n_basap])
input.metadata$population <- factor(input.metadata$population, pop.ord)
input.metadata$Age_bin <- factor(input.metadata$Age_bin, levels = c("18 and above", "14 and below"))
#input.metadata$Country <- factor(input.metadata$Country, levels =c("USA","Tanzania","Nepal","Brazil","Cameroon","Peru","Malaysia"))
input.metadata$Author <- factor(input.metadata$Author, levels =c("Carter","Conteville","Rampelli","Obregon-Tito","Tee"))
input.metadata$Lifestyle <- factor(input.metadata$Lifestyle, levels =c("Industrial", "H-G"))
input.metadata$Geographic_region <- factor(input.metadata$Geographic_region, levels =c("America","Asia"))
input.metadata$Hadza_pop <- factor(input.metadata$Hadza_pop, levels =c("0","1"))
print("California was used as the reference population.")
print("18 and above was used as the reference age.")
write.table(input.comm, "data/maaslin3_input_comm.tsv",quote = F,sep = "\t")
write.table(input.metadata,"data/maaslin3_input_metadata.tsv",quote = F,sep = "\t", row.names = FALSE)


#########################
# fit univariate models #
#########################

# population as fixed effect

output_dir0 <- "Results/Maaslin3/Fixed_random/"
output_dir <- paste0(output_dir0,"HG_pop_count_geo_ran_0")

fit <- maaslin3(input_data = input.comm,
                input_metadata = input.metadata,
                output=output_dir,
                fixed_effects = c( "population","N_reads_clean", "Geographic_region" ),
                normalization = "CLR",
                transform = "None",
                save_models=T,
                warn_prevalence=FALSE,
                small_random_effects= FALSE,
               # random_effects = c("Author" ),
                cores = 3)



input.file <- read.delim("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/all_results.tsv")
input.file <- input.file[input.file$qval_individual<0.05, ]
input.file <- input.file[!is.na(input.file$qval_individual), ]
input.file <- input.file[input.file$metadata== "Lifestyle", ]
signif.features <- input.file



########################Look at species significant 

signficant.species <- merge(signif.features, full.taxonomy, by ="feature")

write.csv(x = signficant.species, "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/significant_species.csv")
#######################
## MAASLIN author ranef ##
#######################
library(lme4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load model 
models <- readRDS("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/fits/models_logistic.rds")

pop.varex <- data.frame()
for (i in names(models)){
  test <- models[[i]]
  ranef.var <- data.frame(VarCorr(test))[,c("grp","var1","vcov","sdcor")]
  varex <- ranef.var$vcov
  names(varex) <- ranef.var$grp
  
  X <- model.matrix(test)
  beta <-fixef(test)
  eta <- X %*% beta
  var_fixed <- var(as.numeric(eta))
  var_residual <- (pi^2)/3
  var_random <- sum(varex)
  var_total <- var_fixed + var_random + var_residual
  varex <- round(varex*100/var_total,1)
  
  tmp <- data.frame(feature=i,pop_var=varex[1],author_var = varex[2],
                    fixed_var=round(100*var_fixed/var_total,1),
                    total_ranef=sum(varex),
                    row.names = NULL)
  
  pop.varex <- rbind(pop.varex,tmp)
}



# most abundant in
#df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund <- input.comm
df.abund$sample_id <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sample_id",variable.name = "feature",value.name = "CPM")

#input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata_2<- input.metadata
#input.metadata_2$sample_id <- rownames(input.metadata_2)
input.metadata_2 <- input.metadata_2[,c("sample_id","population")]
full.taxonomy <- Abund_strong_2 %>%
  select("Genome", "phylum","class","order","family","genus","species") %>%
  rename(feature = Genome)
df.abund <- merge.data.frame(input.metadata_2,df.abund, by="sample_id")

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
            "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/all_results_popvar_author_ranef.tsv", 
            quote = F,sep = "\t", row.names = F)

popvar.cutoff <- 35 # in percentage
strong.popranef <- unique(popranef.df[popranef.df$total_ranef>=popvar.cutoff,"feature"])
length(strong.popranef)


##################
## diet volcano ##
##################
df <- read.delim("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop//all_results.tsv")
#df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_noranef/all_results_wTax.tsv")
length(unique(df$feature)) -> ct
df <- subset(df, qval_individual < 0.05)
length(unique(df$feature)) -> sig.ct

sig.ct *100 / ct
df$Genome <- df$feature
df <- merge(df, Phylum, by= "feature" )

ggdata <- df
ggdata$association <- ggdata$value
ggdata$association[ggdata$qval_individual>=0.05] <- "none"
#ggdata$association <- factor(ggdata$association, c("Orang Asli","California","Nepali","Baka","Matses","Yanomami","Norman","Urban malaysian","none"))
ggdata$association <- factor(ggdata$association, c("H-G"))

write.table(ggdata ,"Results/Maaslin/population_author_ranef/allresults_info.txt",
            quote = F, row.names = F, sep = "\t")
#ggdata <- df[,c("feature","value","coef","N.not.0","pval","qval","species","Genome")]
ggdata <- ggdata[order(abs(ggdata$coef),decreasing = T),]
# ggdata$label <- seq(1,nrow(ggdata),by=1)
ggdata$label <- ggdata$species
ggdata$label[ggdata$qval_individual>=0.05] <- ""
rownames(ggdata) <- NULL
ggdata <- ggdata[ggdata$metadata == "Lifestyle",]

diet.label <- ggdata[ggdata$label!="",]
diet.label <- merge.data.frame(diet.label, full.taxonomy[,c("feature","phylum","family")], by="feature")
rownames(diet.label) <- NULL
write.table(diet.label,"Results/Maaslin/population_author_ranef/all_results_pop_label.tsv",
            quote = F, row.names = F, sep = "\t")

volc <- ggplot(ggdata, aes(y=-log10(qval_individual), x=coef, fill=association)) + theme_bw() +
  geom_text(aes(label=label), col="grey70", size=3,nudge_y =0.1) +
  annotate("text", x = 0, y = -log10(0.045), label = "qval=0.05", vjust=0, col="red")+
  geom_hline(yintercept = -log10(0.05), col="red", lty=2)+
  geom_point(size=5, pch=21, col="grey30") +
  scale_fill_manual(values = c("#2E8B57","#FF69B4","#63B8FF","grey70","#20B2AA","#8B008B","#FF7F50","#DC143C","#FFD700")) +
  scale_x_continuous(limits = c(-6,7), breaks = sort(c(seq(-4,4, by=1),0)))+
  scale_y_continuous(limits=c(1.5, 6)) +
  labs(title = "Lifestyle Correlations (Multi-variate)",
       caption = "q-val is p-val adjusted with Benjamin Hochberg method, q-val<0.05 considered as significant") +
  ylab("Significance -log10(qval)") + xlab("Coefficient")

volc

outpath <- "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/"
for(d in c("png","pdf")){
  ggsave(plot = volc, 
         filename = paste0(outpath,"lifestyle_volcano_2.",d), scale = 1.5,
         device = d,width = 10,height = 6,units = "in",dpi = "print")}


# Compare counts
#significant.count <- readRDS("Results/Maaslin/population_author_ranef/significant_results.tsv")
significant.count <- read.delim("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/significant_results.tsv")
library(ggplot2)
library(ggpubr)
library(ggvenn)

# ############################################################



###############
## diet list significant associations  ##
###############
df <- read.csv("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/all_results.tsv", header=T,sep= "\t")
ggdata<- df
#ggdata <- df[,c("feature","value","coef","pval","N.not.0","qval","phylum","class","family","species")]
ggdata <- ggdata[ggdata$qval_individual<0.05,]
ggdata <- ggdata[ggdata$metadata == "Lifestyle",]
ggdata <-  ggdata[!is.na(ggdata$feature), ]
diet.MAGs <- ggdata$feature
length(diet.MAGs)

# order data
ggdata <- ggdata[order(ggdata$coef,decreasing = F),]
#ggdata$value <- factor(ggdata$value, c("Orang Asli","California","Nepali","Baka","Matses","Yanomami","Norman","Urban malaysian","none"))
ggdata$value <- factor(ggdata$value, c("H-G"))

mag.levels <- as.character(unique(ggdata$feature))
ggdata$feature <- factor(ggdata$feature, mag.levels)

#phylum.tree.order <- rev(readRDS("output_files/revised_MAG/phylum.tree.order.rds"))
phylum.tree.order <-unique(full.taxonomy$phylum)
ggdata <- merge(ggdata, full.taxonomy, by="feature")
#phylum.tree.order <- c(phylum.tree.order,"Methanobacteriota")
ggdata$phylum <- factor(ggdata$phylum, phylum.tree.order)

diet.list.df <- ggdata 

########
# plot #
########


library(reshape2)
library(dplyr)
#df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund <- input.comm
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

#input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
#input.metadata_2$sample_id <- rownames(input.metadata)
input.metadata_2<- input.metadata
input.metadata_2 <- input.metadata[,c("sample_id","population")]
input.metadata_2 <- input.metadata_2 %>%
  rename(sampleid = sample_id)
df.abund <- merge.data.frame(input.metadata_2,df.abund, by="sampleid")
df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

n <- length(unique(df.abund$sampleid))


df.abund %>% group_by(feature, phylum) %>% summarise(abundance=sum(CPM)/(n*10^6)) -> total_abund

heat.df.total <- total_abund[total_abund$feature %in% mag.levels,]
heat.df.total$feature <- factor(heat.df.total$feature, mag.levels)
heat.df.total$phylum <- factor(heat.df.total$phylum, phylum.tree.order)


df.abund %>% group_by(population,feature, phylum) %>% 
  summarise(n=length(unique(sampleid)),
            popCPM=sum(CPM),
            popMedian=median(CPM),
            popMean=mean(CPM),
            popSD=sd(CPM)) -> pop_abund



df.abund %>% group_by(feature) %>% summarise(totalCPM=sum(CPM)) -> total_abund

heat.df.pop <- merge.data.frame(pop_abund, total_abund, by="feature")

pop_abund %>% group_by(feature) %>% summarise(sum.popMean=sum(popMean)) -> sumpopMean
heat.df.pop <- merge.data.frame(heat.df.pop, sumpopMean, by="feature")

heat.df.pop <- heat.df.pop[order(heat.df.pop$feature),]


heat.df.pop$normalised.abund <- 100*heat.df.pop$popCPM/heat.df.pop$totalCPM

heat.df.pop$normalised.average <- 100*heat.df.pop$popMean/heat.df.pop$sum.popMean

# order data
heat.df.pop <- heat.df.pop[heat.df.pop$feature %in% mag.levels,]
heat.df.pop$feature <- factor(heat.df.pop$feature,mag.levels)
heat.df.pop$phylum <- factor(heat.df.pop$phylum, phylum.tree.order)
heat.df.pop$population <- factor(heat.df.pop$population,
                                 c("Hadza","Orang Asli","California","Nepali","Baka",
                                   "Matses","Yanomami","Norman","Urban malaysian"))


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
        text=element_text(size=20), axis.text=element_text(size=20),
        #axis.text.y = element_text(angle=0, hjust = 1, size=8),
        axis.text.x=element_text(colour = "black", size=16),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y= element_blank())+ 
  scale_color_manual(values = "pink") +
  scale_x_continuous(breaks = seq(-10,10, by=1),limits = c(-5,7), expand = c(0.06,0.06))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Coefficients") +
  labs(title="MAGs with Significant Lifestyle Associations (qval<0.05)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_point(aes(col= "pink"),size=3) +
  geom_text(aes(label=species), angle=0, nudge_x = 0.1, 
            size=3, col="grey40", hjust=0, vjust=0.2)
ggsave("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/MAGs_plot.png", a, width = 10, height = 20, dpi = 300)

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
        axis.text.x = element_text(size=16, colour = "black",angle=60, hjust = 1),
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
        axis.text.x = element_text(size=16, colour = "black",angle=60, hjust = 1),
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
ggsave("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/MAGs_plot.png", a, width = 10, height = 20, dpi = 300)

ggsave(filename = "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/MAGs_final_plot.pdf",diet.list,device = "pdf",width = 15, height = 20,units = "in",dpi = 300 )
saveRDS(grid.arrange(a,b,c,layout_matrix=lay),
        file = "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/MAGs_final_plot.rds")



############################### plot variance 
### POP RANEF ####
# create figure data frame
ggdata <- popranef.df
ggdata$topPop <- factor(ggdata$topPop, c("California","Baka","Hadza","Nepali", "Norman","Matses","Urban malaysian", "Orang Asli", "Yanomami"))
ggdata <- ggdata[ggdata$total_ranef!=0,]
ggdata$label <- ggdata$feature

ggdata <- merge.data.frame(ggdata, full.taxonomy, by="feature")
write.table(ggdata,"Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/model_popvar.tsv", quote = F, row.names = F, sep = "\t")


summary(ggdata$total_ranef)
summary(ggdata$author_var)
quantile(ggdata$total_ranef)
quantile(ggdata$author_var)
hist(ggdata$total_ranef, breaks=30)
hist(ggdata$author_var, breaks=30)

ggdata %>% summarise(n=length(total_ranef), 
                     average=mean(total_ranef),
                     sd=sd(total_ranef),
                     se=sd(total_ranef)/length(total_ranef),
                     min=min(total_ranef),
                     Q1=quantile(total_ranef)[2],
                     med=median(total_ranef),
                     Q3=quantile(total_ranef)[4],
                     max=max(total_ranef))

ggdata %>% 
  filter(total_ranef > 35) %>% 
  group_by(phylum) %>%
  summarise(n=length(total_ranef), 
            average=mean(total_ranef),
            sd=sd(total_ranef),
            se=sd(total_ranef)/length(total_ranef),
            min=min(total_ranef),
            Q1=quantile(total_ranef)[2],
            med=median(total_ranef),
            Q3=quantile(total_ranef)[4],
            max=max(total_ranef))



# Filter MAGs 
ggdata <- ggdata[order(ggdata$total_ranef, decreasing = T),]
best.pop <- ggdata[ggdata$total_ranef>35,] # varex > 35%
# best.pop <- ggdata 

# set phylum colours
#phylum.cols <- read.delim("input_files/phylum_colours_bar.tsv", header=T)
phylum.cols <- c(
  "#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#ffbb78",
  "#98df8a", "#aec7e8", "#ff9896", "#9467bd", "#c5b0d5",
  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#8B008B",
  "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
)
names(phylum.cols) <- c("Bacteroidota","Bacillota_A","Spirochaetota","Actinomycetota","Pseudomonadota","Verrucomicrobiota","Bacillota","Bacillota_B","Bacillota_C","Cyanobacteriota","Elusimicrobiota","Methanobacteriota","Desulfobacterota","Campylobacterota","Fibrobacterota","Myxococcota","Thermoplasmatota","Eremiobacterota","Planctomycetota","Fusobacteriota")

# tmp <- phylum.cols$cols
# names(tmp) <- phylum.cols$phylum
# phylum.cols <- tmp

# # set population colours
# popcols <- readRDS("input_files/pop_mycols.rds")
# popcols <- popcols[names(popcols) %in% unique(input.metadata$population)]
# 
# table(best.pop$phylum)

# order data
best.pop <- best.pop[order(best.pop$total_ranef,decreasing = T),]
best.pop$feature <- factor(best.pop$feature, level=unique(best.pop$feature))
mag.levels <- unique(best.pop$feature)
best.pop$feature <- factor(best.pop$feature, rev(mag.levels))
best.pop <- best.pop[order(rev(mag.levels)),]

# phylum.tree.order <- rev(readRDS("output_files/revised_MAG/phylum.tree.order.rds"))
# phylum.tree.order <- c(phylum.tree.order,"Methanobacteriota")
best.pop$phylum <- factor(best.pop$phylum, c("Bacteroidota","Bacillota_A","Spirochaetota","Actinomycetota","Pseudomonadota","Verrucomicrobiota","Bacillota","Bacillota_B","Bacillota_C","Cyanobacteriota","Elusimicrobiota","Methanobacteriota","Desulfobacterota","Campylobacterota","Fibrobacterota","Myxococcota","Thermoplasmatota","Eremiobacterota","Planctomycetota","Fusobacteriota"))

rownames(best.pop) <- NULL

# add diet associations
df <- read.csv("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/significant_species.csv", sep= ";")
diet.assoc <- df[df$qval_individual<0.05, c("feature","value")]
best.pop$diet.assoc <- "population"
for(i in 1:nrow(best.pop)){
  x <- as.character(best.pop$feature)[i]
  y <- x %in% diet.assoc$feature
  z <- diet.assoc$value[diet.assoc$feature==x]
  if(isTRUE(y)){
    best.pop[i,"diet.assoc"] <- y
  }
}

########
# plot #
########
diet.cols <- c("#00A600","#FF69B4","#0000FF","#7A67EE","grey70")
names(diet.cols) <- c("Sago","Pork","Chicken","population","not significant")

#ggplot(best.pop, aes(x=pop_var,y=feature))+ theme_classic() +
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.background = element_rect(color = "black",fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.text.y = element_text(angle=0, hjust = 1, size=8),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.title.y= element_blank(),
#        panel.grid.major.x = element_line(colour="grey",linetype = 2)) + 
#  scale_color_manual(values = diet.cols) +
#  scale_x_continuous(breaks = seq(0,100, by=10),limits = c(35,105), expand = c(0,0),labels = paste0(seq(0,100, by=10),"%"))+
#  scale_y_discrete(expand = c(0,0.5))+
#  xlab("Variance Explained") +
#  labs(title="MAGs with Strongest Population Effect (varience explained >35%)")+
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_point(aes(col=diet.assoc),size=3) +
#  geom_text(aes(label=species), angle=0, nudge_x = 0.5, 
#            size=3, col="grey40", hjust=0, vjust=0.2)

# Total Abundance
mag.levels <- as.character(unique(best.pop$feature))

df.abund <- read.delim("data/maaslin3_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

# input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
# input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sample_id","population")]
input.metadata <- input.metadata %>%
  rename(sampleid =sample_id)
df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")

df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(feature, phylum) %>% summarise(abundance=sum(CPM)/(n*10^6)) -> total_abund

# order data
heat.df.total <- total_abund[total_abund$feature %in% as.character(mag.levels),]
heat.df.total$feature <- factor(heat.df.total$feature, mag.levels)
heat.df.total$phylum <- factor(heat.df.total$phylum, phylum.tree.order)

########
# plot #
########

# ggplot(heat.df.total)+ theme_classic() +
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_tile(aes(x="%Abundance (log10)", y=feature,fill=log(abundance)))+
#  scale_fill_gradient(low = "grey90", high = "black",)+
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.ticks.y= element_blank(),
#        axis.title.y= element_blank()) +
#  labs(title="Total Abundance") +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0.5))+
#   xlab("")


# add heatmap abundance # normalised average
df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")

df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(population,feature, phylum) %>% 
  summarise(popCPM=sum(CPM),
            popMedian=median(CPM),
            popMean=mean(CPM),
            popSD=sd(CPM)) -> pop_abund



df.abund %>% group_by(feature) %>% summarise(totalCPM=sum(CPM)) -> total_abund

heat.df.pop <- merge.data.frame(pop_abund, total_abund, by="feature")

pop_abund %>% group_by(feature) %>% summarise(sum.popMean=sum(popMean)) -> sumpopMean
heat.df.pop <- merge.data.frame(heat.df.pop, sumpopMean, by="feature")

heat.df.pop <- heat.df.pop[order(heat.df.pop$feature),]

heat.df.pop$normalised.abund <- 100*heat.df.pop$popCPM/heat.df.pop$totalCPM
heat.df.pop$normalised.average <- 100*heat.df.pop$popMean/heat.df.pop$sum.popMean

# order data
heat.df.pop <- heat.df.pop[heat.df.pop$feature %in% as.character(mag.levels),]
heat.df.pop$feature <- factor(heat.df.pop$feature,mag.levels)
heat.df.pop$phylum <- factor(heat.df.pop$phylum, phylum.tree.order)
heat.df.pop$population <- factor(heat.df.pop$population,
                                 c("Hadza","Baka","Matses","Yanomami","Nepali","Orang Asli","Urban malaysian","Norman","California"))

########
# plot #
########

# ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  theme_minimal() + 
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.ticks.y= element_blank(),
#        axis.title.y= element_blank()) +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0.5))+
#  geom_tile()+ scale_fill_gradient(low="white",high="black") +
#  labs(title="Normalised Average") +
#  xlab("")

# combine graph
a <- ggplot(best.pop, aes(x=total_ranef,y=feature,unit = "px"))+ 
  theme_minimal() + 
  theme(plot.margin = margin(0,0,17,0, unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        panel.grid.minor.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_rect(fill="white",linewidth = 0),
        text=element_text(size=20), axis.text=element_text(size=20),
        #axis.text.y = element_text(angle=0, hjust = 1, size=8),
        axis.text.x=element_text(colour = "black", size=16),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y= element_blank())+
  scale_color_manual(values = diet.cols) +
  scale_x_continuous(breaks = seq(0,100, by=10),limits = c(35,105), 
                     expand = c(0.06,0.06),labels = paste0(seq(0,100, by=10),"%"))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Variance Explained") +
  labs(title="MAGs with Population, Author Effect (variance explained >35%)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_point(aes(col="black"),size=3)+
  geom_text(aes(label=species), angle=0, nudge_x = 1.2, 
            size=5, col="grey40", hjust=0, vjust=0.5)

a
ggsave(filename = "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/MAGs_variance_plot.pdf",a,device = "pdf",width = 15, height = 20,units = "in",dpi = 300 )
ggsave("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/MAGs_variance_plot.pdf", a, width = 10, height = 20, dpi = 300)
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
        panel.border = element_rect(colour = "black", fill=NA),
        # strip.text.y.left = element_text(angle=0),
        # strip.background = element_rect(fill="white",linewidth = 0),
        strip.text = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=20), axis.text=element_text(size=20),
        axis.ticks.x= element_blank(),
        axis.text.x=element_text(colour = "black", size=16, angle = 60, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        #axis.text.y = element_text(size=8), axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  labs(title="") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("")
b
c <- ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  theme_minimal() + 
  theme(plot.margin = margin(0,0,3,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        # strip.text.y.left = element_text(angle=0),
        # strip.background = element_rect(fill="white",linewidth = 0),
        strip.text = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=20), axis.text=element_text(size=20),
        axis.ticks.x= element_blank(),
        axis.text.x=element_text(colour = "black", size=16,angle = 60, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  geom_tile()+ scale_fill_gradient(low="grey95",high="black") +
  labs(title="") +
  xlab("")
c
lay <- rbind(c(rep(1,33),
               rep(2,1),
               rep(3,4)))
# draw plot
popvar.list <- grid.arrange(a,b,c,layout_matrix=lay)
ggsave("Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/variances_plot.png", a, width = 10, height = 20, dpi = 300)
ggsave(filename = "Results/Maaslin3/Fixed_random/Lifestyle_count__geo_ran_author_pop/figures/variances_final_plot.pdf",popvar.list,device = "pdf",width = 15, height = 20,units = "in",dpi = 300 )
saveRDS(popvar.list, "output_files/revised_MAG/Figure3/Fig3F.rds")













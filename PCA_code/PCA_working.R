# This script performs and plots a PCA of microbiome variation based relative abundances. 

# Load libraries
library(vegan)
library(compositions)

# Set working directory
setwd("/Users/hl636/Documents/Hirzi/Cambridge/Student projects/Lucie/summary/")

# Define ordination method: "clr_pca" or "aitchison_capscale". The former is written manually, whereas the latter calls functions vegdist and capscale from the package "vegan".
# Note that clr_pca (i.e., performing an unscaled PCA on CLR-transformed compositional data (i.e. relative abundance table)) is identical to applying PCoA on the Aitchison distance matrix (i.e., "aitchison_capscale").
ord_method <- "clr_pca" # "aitchison_capscale" or "clr_pca"
if(ord_method == "clr_pca") {
  scale_PCA <- FALSE
} else if(ord_method == "aitchison_capscale") {
  distance_metric <- "robust.aitchison" # "jaccard", "aitchison", "robust.aitchison" or any vegdist method. "aitchison" and "robust.aitchison" give identical results.
}

# Read in data
#mypop <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Student projects/Lucie/ModernHG_popIndo.txt", header=TRUE)
#mypop <- sample_pop_all_list
mypop <- read.csv("Sample_pop_meta/metadata_pop_all_final.csv", header= TRUE, sep = ";")
samples_exclude <- c("ERR7738433", "ERR7738234", "ERR7738665") # ERR7738433_1.fastq.gz and ERR7738234_2.fastq.gz were not downloaded successfully, hence all downstream steps with this sample failed. ERR7738665 paired reads have different names.
sampleid <- mypop$sample_id [! mypop$sample_id %in% samples_exclude]

df <- read.delim("data/bwa_counts_total_filtered_wMetadata_15mil_s95_strongFilter_relAbund.tsv")

# Filter by MAGs (in case you e.g. want to perform for rare MAGs)
df_filt <- df

# Filter by sample
comm.df <- df_filt[,colnames(df_filt) %in% sampleid]
#idx_firstSample <- match("centrality", colnames(df)) + 1
comm.df <- comm.df*10^6
rownames(comm.df) <- df_filt$Genome
comm.df <- data.frame(t(comm.df))
comm.df <- comm.df[sampleid,]
comm.df <- comm.df[rowSums(comm.df)>0,]
table(colSums(comm.df)==0)
dim(comm.df)
comm.df<- na.omit(comm.df)
dim(comm.df)
pc.sampleid <- rownames(comm.df)

# metadata
env.df <- mypop[mypop$sample_id %in% sampleid,]
removed.samples <- mypop$sample_id[!mypop$sample_id %in% sampleid]
rownames(env.df) <- env.df$sample_id
env.df <- env.df[,c("sample_id","population")]
env.df <- env.df[sampleid,]
data.frame(table(env.df$population))

############
### PCA ###
############

if(ord_method == "clr_pca") {
  rowSums(comm.df) # should be ~1x10ˆ6
  # log(0) is undefined (= -Inf). Note that compositions::clr return 0 for log(0) through if/else statement, however, we want to avoid this behavior (since it is incorrect). 
  # To work around, we add an offset to the data (here, an order of magnitude smaller than the dataset's minimum value), or 1 for counts data. For reference, see: http://mixomics.org/mixmc/mixmc-preprocessing/
  comm.df.clr <- as.data.frame(t(clr(comm.df + 1)))
  comm.df.clr.t <- as.data.frame(t(comm.df.clr))
  if(scale_PCA == TRUE) {
    pc <- prcomp(comm.df.clr.t, center = TRUE, scale. = TRUE)
  } else if (scale_PCA == FALSE) {
    pc <- prcomp(comm.df.clr.t, center = TRUE, scale. = FALSE)
  }
  pc_summary <- summary(pc)
  pc_df <- as.data.frame(pc$x)
  var <- pc_summary$importance[2,]
} else if(ord_method == "aitchison_capscale") {
  comm.dist <- vegdist(comm.df+1, method=distance_metric)
  cap <- capscale(comm.dist~1, comm=comm.df)
  eig <- cap$CA$eig
  var <- eig/sum(eig)
}

# presence absence filter
# raup.df <- comm.df
# raup.df <- raup.df[rowSums(raup.df) > 0,]
# raup.df[raup.df>0] <- 1
# comm.dist <- vegdist(raup.df, method="raup")
# pc.sampleid <- rownames(raup.df)
# dim(raup.df)

####################
# Draw pretty plot #
####################

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(randomcoloR)

# set colours
mycols <- data.frame(mycols = distinctColorPalette(length(unique(mypop$population))), group = unique(mypop$population))

# select PC
PC <- c(1,2)
###########
# prepare plot data for collation
if(ord_method == "clr_pca") {
  xycoord <- cbind(pc_df[, c(1,2)], rownames(pc_df))
  axis_labs <- paste("PC",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
} else if(ord_method == "aitchison_capscale") {
  xycoord <- summary(cap)
  xycoord <- data.frame(scores(cap, display = "sites"))
  #xycoord <- data.frame(xycoord$sites[,PC])
  xycoord$sample_id <- rownames(xycoord)
  axis_labs <- paste("PCo",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
}
rownames(xycoord) <- NULL
pco.coord <- xycoord # save for later
colnames(xycoord) <- c("Xaxis","Yaxis","sample_id")

# limit settings normal
#lims <- c(-6,4); ylims <- c(-5,5)
xlims <- c(min(xycoord$Xaxis)*1.2,max(xycoord$Xaxis)*1.2); ylims <- c(min(xycoord$Yaxis)*1.2,max(xycoord$Yaxis)*1.2)
# xylimits settings presence/absence
# xlims <- c(-2,2); ylims <- c(-2,2)

# set population order
pop.ord <- mycols$group

# set environmental variables
ggdata <- merge.data.frame(xycoord, mypop, by="sample_id")

ggdata$population <- factor(ggdata$population, pop.ord)
ggdata <- ggdata[order(ggdata$population),]

ggdata$sampleid <- factor(ggdata$sample_id, unique(ggdata$sample_id))
rownames(ggdata) <-  NULL

# set axis title
xtitle <- axis_labs[PC[1]]
ytitle <- axis_labs[PC[2]]

# set color by country
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord
mycols <- c(mycols)

# make base plot
ggdata$population <- factor(ggdata$population, pop.ord)
a <- ggplot() +
  theme_bw() +
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept=0, lty="dashed")

# draw the rest of the plot
a <- a +
  geom_point(data=ggdata, aes(x=Xaxis, y=Yaxis, fill=population), col="black", size=3, pch=21) +
  geom_point(data=ggdata, aes(x=Xaxis, y=Yaxis, fill=population), col="black", size=3, pch=21) +
  scale_fill_manual(values=mycols) +
  scale_x_continuous(expand = c(0,0),
                     limits = xlims,
                     breaks = round(seq(xlims[1],xlims[2],length.out = 5)),
                     position = "top")+
  scale_y_continuous(expand = c(0,0),
                     limits = ylims,
                     breaks = round(seq(ylims[1],ylims[2],length.out = 5)),
                     position = "left") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        plot.margin=margin(0,0,0,0,"cm"),
        legend.position = "none")+
  xlab(xtitle) + ylab(ytitle)


# PCo Group Boxplot X axis
pop_cols <- mycols[pop.ord] # set color by population

ggdata$population <- factor(ggdata$population, rev(pop.ord))
ggdata <- ggdata[order(ggdata$population),]
b <- ggplot(ggdata,
            aes(x=population, y=Xaxis, fill=population))+
  theme_bw()+
  geom_hline(yintercept = 0, lty="dashed") +
  geom_boxplot(width=0.8, outlier.size = 2, outlier.shape=21,
               outlier.fill = "white") +
  scale_fill_manual(values=pop_cols[pop.ord]) +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_blank(),
        axis.title.y = element_text(colour = "white",size=12),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin=margin(0,0,0,0,"cm"),
        legend.position="none")+
  scale_y_continuous(expand = c(0,0),
                     limits = xlims,
                     breaks = round(seq(xlims[1],xlims[2],length.out = 5)),
                     position = "left") +
  xlab("sample distance") + ylab(" ")+
  coord_flip()

# PCo Group Boxplot Y axis
ggdata$population <- factor(ggdata$population, pop.ord)
c <- ggplot(ggdata,
            aes(x=population, y=Yaxis, fill=population)) +
  theme_bw()+
  geom_hline(yintercept = 0, lty="dashed")+
  geom_boxplot(width=0.8, outlier.size = 2,
               outlier.shape=21, outlier.fill = "white") +
  scale_fill_manual(values=pop_cols) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_text(colour = "white",size=12),
        plot.margin=margin(0,0,0,0,"cm"),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),
                     limits = ylims,
                     breaks = round(seq(ylims[1],ylims[2],length.out = 5)),
                     position = "right")+
  scale_x_discrete(position="top")+
  xlab("sample distance") + ylab(" ")


# fix labels and extract legend

pop_count <- data.frame(table(mypop$population))
colnames(pop_count)  <- c("population","sample")
pop_count$population <- factor(pop_count$population, pop.ord)

poplab <-paste0(pop_count$population," (n=",pop_count$sample,")")
names(poplab) <- pop_count$population

ggdata$label <- c()
for(pop in names(poplab)){
  ggdata[ggdata$population == pop, "label"] <- poplab[pop]
}
ggdata <- ggdata[order(ggdata$population),]
ggdata$label <- factor(ggdata$label, unique(ggdata$label))

legend.col <-  mycols[as.character(unique(ggdata$population))]
names(legend.col) <- poplab[names(legend.col)]

p <- ggplot(ggdata, aes(x=label, y=Xaxis, fill=label))+ theme_bw() +
  geom_boxplot(width=0.8, outlier.size = 2, outlier.shape=21, outlier.fill = "white") +
  scale_fill_manual(values=legend.col) +
  theme(legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.background = element_blank())
legend <- cowplot::get_legend(p)
legend <- as_ggplot(legend) + theme(plot.margin = margin(-1.5,1,0,0,"cm"))


# set plot layout
lay <- rbind(c(1,1,1,1,2,2),
             c(1,1,1,1,2,2),
             c(1,1,1,1,2,2),
             c(1,1,1,1,2,2),
             c(3,3,3,3,4,4),
             c(3,3,3,3,4,4))
# draw plot
print("Plot: Aitchison Distance for CLR-transformed Data")
grid.arrange(a,c,b,legend, layout_matrix=lay)



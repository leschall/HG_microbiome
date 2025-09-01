###Venn diagram``
library(ggvenn)

## Venn Analysis ###
diet.cols <- c("#2E8B57","#DB7093","#63B8FF","#DC143C")
#names(diet.cols) <- c("Sago","Pork","Chicken")

multivar.df <- read.csv("Results/Maaslin3/Fixed/significant_results_all.csv", sep = ";" )

vs <- list(Lifestyle=multivar.df$feature[multivar.df$metadata=="Lifestyle"],
           Read_count=multivar.df$feature[multivar.df$metadata=="N_reads_clean"],
           Continent=multivar.df$feature[multivar.df$metadata=="Geographic_region"],
           Population=multivar.df$feature[multivar.df$metadata=="population"])

p <- ggvenn(vs,text_size = 4, set_name_size = 5, fill_color = diet.cols) +
  coord_fixed() +                      # <— garde des cercles parfaits
  theme(aspect.ratio = 1)              # (optionnel) renforce le ratio carré
print(p)

ggsave(p, file="Figures/final/venn_fixed.png",device = "png",width = 8,height = 8)
ggsave(p, file="Figures/final/venn_fixed.pdf",device = "pdf",width = 8,height = 8)
saveRDS(p, file="output_files/revised_MAG/Figure3/Fig3C.rds")



vs2 <- list(Population=MAG.list[["population_noranef"]],
            Diet=MAG.list[["Sago+Chicken+Pork_noranef"]])

p2 <- ggvenn(vs2,text_size = 4,set_name_size = 5, fill_color = c("#7A67EE","#FA8072"))
print(p2)

ggsave(p2, file="output_files/revised_MAG/Figure3/Fig3D.svg",device = "svg",width = 8,height = 8)
ggsave(p2, file="output_files/revised_MAG/Figure3/Fig3D.pdf",device = "pdf",width = 8,height = 8)
saveRDS(p2, file="output_files/revised_MAG/Figure3/Fig3D.rds")


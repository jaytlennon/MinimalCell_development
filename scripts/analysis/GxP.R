rm(list = ls())
#library(ggplot2);library(tidyverse);library(dplyr);library(ggpubr);library(Hmisc)
package.list <- c('vegan', 'ade4', 'viridis', 'gplots', 'BiodiversityR', 'indicspecies', 'ggplot2', 'tidyverse', 'dplyr', 'ggpubr', 'Hmisc')
for (package in package.list) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


myt <- read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\GxP_Mm_essential.csv")
myt_ALL.shared<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\GxP_Mm_all.shared.csv")

myt <- read.csv("data2/GxP_Mm_essential.csv")
myt_ALL.shared <- read.csv("data2/GxP_Mm_all.shared.csv")


myt$sample <- as.factor(myt$sample)
myt$treatment <- as.factor(myt$treatment)
levels(myt$treatment)[levels(myt$treatment)=="s1"] <- "Wildtype"
levels(myt$treatment)[levels(myt$treatment)=="s3"] <- "Minimal"




#data(doubs)
#typeof(doubs)
#typeof(doubs$fish)

mygxp<-as_tibble(myt)
gxpnum<-mygxp[-1]
gxpnum<-gxpnum[-1]
gxpnum<-data.matrix(gxpnum)

adonis(gxpnum ~ myt$treatment, method = "bray", permutations = 99999)
#same result with Jaccard, Soerensen, Horn, Kulczynski, &c

gxp_BC_triangle<-vegdist(gxpnum,method = 'bray')
gxp_BC_double<-vegdist(gxpnum,method='bray',upper=TRUE,diag=TRUE)#again the results are similar when using other distance metrics
gxp_BC<-as_tibble(data.matrix(gxp_BC_double))

gxp_BC_meta<-add_column(gxp_BC, sample = myt$sample, .before = 1)
gxp_BC_meta<-add_column(gxp_BC_meta, treatment = myt$treatment, .before = 2)
gxp_BC



gxp.pcoa.eig <- cmdscale(gxp_BC, eig = TRUE, k = 3)
explainvar1 <- round(gxp.pcoa.eig$eig[1] / sum(gxp.pcoa.eig$eig), 4) * 100#4 is specifying the number of decimal places
explainvar2 <- round(gxp.pcoa.eig$eig[2] / sum(gxp.pcoa.eig$eig), 4) * 100
explainvar3 <- round(gxp.pcoa.eig$eig[3] / sum(gxp.pcoa.eig$eig), 4) * 100
# so here simply the Euclidean distance is used, see distances above
print(explainvar1)#42.12%
pco1lab<-toString(explainvar1)
print(explainvar2)#16.46%
pco2lab<-toString(explainvar2)
print(explainvar3)#15.14%


gxp_pcoa<- cmdscale(gxp_BC, k =3)
gxp.pcoa <- as.data.frame(gxp_pcoa)
gxp.pcoa$treatment <- myt$treatment
row.names(gxp.pcoa)<-myt$sample
gxp.pcoa$sample<-myt$sample
names(gxp.pcoa)[1:3] <- c('PC1', 'PC2', 'PC3')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "blue", "red",
               "#0072B2", "#D55E00", "#CC79A7")
gxpplot <- ggplot(gxp.pcoa, aes(x = PC1, y = PC3, colour = treatment, 
                                 label = row.names(gxp.pcoa)))
gxpplot + geom_point(size =5) +
  scale_colour_manual(values = cbPalette[4:5]) +
  geom_hline(yintercept = 0, linetype=117, size = 1.5, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype=117, size = 1.5, color = "darkgrey") +
  #geom_text(col = 'black', size = 0)+
  stat_ellipse(level=0.95, show.legend = FALSE)+
  labs(x="\nPCo1 (42.12%)",y="PCo2 (16.46%)\n",col="Genotype") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_line(color="black", size = 3.5)) + 
  theme(legend.key=element_blank(),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x = element_line(color="black", size = 1.5), axis.line.y = element_line(color="black", size = 1.5), axis.ticks.y = element_line(color = "black", size = 1.5), axis.ticks.x = element_blank()) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())

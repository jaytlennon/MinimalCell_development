rm(list = ls())
package.list <- c('ggplot2','tidyverse','dplyr','ggpubr','Hmisc')
for (package in package.list) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
theme_set(
  theme_bw()
)

mydf<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\microscopy\\microbeJ_results\\cases_microbeJ.csv")
#mydf<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\cases_Mm.300_1.2.csv")
#mydf$Strain <- factor(mydf$Strain, levels=c("Wildtype","Minimal"))
mycols<-c("black","black","black","black","black","black","black","black","black","black")
#mycols<-c("white","white","white","white","white","white","white","white","white","white")

areaV<-aov(mydf$Area ~ mydf$Line_exact)
summary(areaV)
TukeyHSD(areaV)


dnds <- ggplot(mydf, aes(x=Line_exact, y=Area))
dnds + geom_jitter(
  aes(shape = Line_exact, color = Line_exact), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = .1, dodge.width = .5),
  size = 1, stroke = 1.3) +
  scale_shape_manual(values = c(1,1,1,1,1,1,1,1,1,1)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = Line_exact),
    fun.data = "mean_se", fun.args = list(mult = (2)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.5, shape=95,
    position = position_dodge(0.2),
    show.legend = FALSE
  ) +
  scale_color_manual(values = mycols) +
  labs(x="\nStrain",y="pixels\n") +
  geom_hline(yintercept = 1, linetype=117, size = 1.25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",legend.key=element_blank(),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line = element_line(colour = "black"), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())
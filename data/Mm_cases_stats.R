rm(list = ls())
#library(ggplot2);library(tidyverse);library(dplyr);library(ggpubr);library(Hmisc)
package.list <- c('vegan', 'ade4', 'viridis', 'gplots', 'BiodiversityR', 'indicspecies', 'ggplot2', 'tidyverse', 'dplyr', 'ggpubr', 'Hmisc')
for (package in package.list) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
idff<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\mean.SE_Mm.300.csv")
BANC<-idff$W_anc1.0[2] / idff$W_anc.own[2]


myt <- read.csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\cases_Mm.300.csv",header=T)
#myt <- read.csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\cases_Mm.300_1.2.csv",header=T)

#
myd<-as_tibble(myt)
my1<-filter(myd,Strain=='syn1.0')
myB<-filter(myd,Strain=='syn3B')
#myd<-add_column(myd, mu=(1/e0))

#A.aov<-aov(myd$A ~ myd$Strain_name)
#summary(A.aov)
#TukeyHSD(A.aov)

V<-var.test(my1$W_anc1.0,myB$W_anc1.0)
V

absoluteW<-t.test(my1$W_anc1.0,myB$W_anc1.0,alternative='two.sided',paired=F,var.equal=T,mu=0)
absoluteW

vi<-var.test(my1$W_anc.own,myB$W_anc.own)
vi

increaseWvsown<-t.test(my1$W_anc.own,myB$W_anc.own,alternative='two.sided',paired=F,var.equal=T,mu=(0))
print(increaseWvsown)

vabsincrease<-var.test(my1$W_anc1.0,myB$W_anc1.0)
vabsincrease


Wabsincrease<-t.test(my1$W_anc1.0,myB$W_anc1.0,alternative="g",paired=F,var.equal=T,mu=(1-BANC))
Wabsincrease
1.3235616-0.8428967


mygg<-ggplot(idff, aes(x=Strain,y=W_anc1.0,fill=Strain)) +
  geom_bar(aes(x=Strain, y=W_anc1.0),position=position_dodge(.1), inherit.aes = TRUE, stat= "identity", colour="black", width = .5) +
  geom_errorbar(aes(Strain, ymin = (W_anc1.0 - 1*W_anc1.0_SE), ymax = (W_anc1.0 + 1*W_anc1.0_SE)),position = position_dodge(1), width = .15, size = 1.5) +
  labs(x="Strain", y="Relative fitness (W)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=22),axis.title=element_text(size=34), axis.line.x = element_line(color="black", size = 1.5), axis.line.y = element_line(color="black", size = 1.5), axis.ticks.y = element_line(color = "black", size = 1.5), axis.ticks.x = element_blank(), axis.ticks.length = unit(0.2, "cm"), aspect.ratio = 2/1) +
  scale_y_continuous(limits = c(0,1.506), expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype=117, size = 1.25) + 
  geom_hline(yintercept = (BANC), linetype=117,col='red', size = 1.25) +
  scale_fill_manual("",values = c("syn1.0"=" dark grey","syn3B"="salmon"))
plot(mygg)

myt
levels(myt$Strain)[levels(NSEOD$Strain)=="3b"] <- "JCVI-syn3B"
levels(NSEOD$Strain)[levels(NSEOD$Strain)=="s1"] <- "JCVI-syn1.0"
levels(NSEOD$Time)[levels(NSEOD$Time)=="anc"] <- "Ancestor"
levels(NSEOD$Time)[levels(NSEOD$Time)=="evolved"] <- "Evolved"
NSEOD$Strain<-factor(NSEOD$Strain, levels=c("JCVI-syn1.0","JCVI-syn3B"))
pgcols<-c("blue","red")
ODfig <- ggplot(NSEevo, aes(x=Strain, y=OD))
ODfig +   geom_errorbar(aes(ymax=c(0.079333333,0.079333333,0.079333333,0.079333333,0.876333333,0.876333333,0.876333333,0.876333333),ymin=c(0.079333333,0.079333333,0.079333333,0.079333333,0.876333333,0.876333333,0.876333333,0.876333333)), color=c("red","red","red","red","blue","blue","blue","blue"),lwd = 2.2,linetype=117) +
  geom_jitter(
    aes(shape = Strain, color = Strain), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
    position = position_jitterdodge(jitter.width = 0.475, dodge.width = 0.2),
    size = 8, stroke=2.3) +
  scale_shape_manual(values = c(0,0)) +
  stat_summary(
    aes(color = Strain),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.5, shape=95,
    position = position_dodge(0.2),
    show.legend=FALSE
  ) +
  scale_color_manual(values = pgcols) +
  labs(x="\nStrain",y="OD 600 after 8 days\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())#+
#theme(legend.position="none")



###########
#Analyze changes in OD during the NSE
myOD<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\OD_600\\20200327_OD.600_summ2way.csv")
NSEOD<-as_tibble(myOD)
NSEOD$Strain <- as.factor(NSEOD$Strain);NSEOD$Time <- as.factor(NSEOD$Time)
levels(NSEOD$Strain)[levels(NSEOD$Strain)=="3b"] <- "JCVI-syn3B"
levels(NSEOD$Strain)[levels(NSEOD$Strain)=="s1"] <- "JCVI-syn1.0"
levels(NSEOD$Time)[levels(NSEOD$Time)=="anc"] <- "Ancestor"
levels(NSEOD$Time)[levels(NSEOD$Time)=="evolved"] <- "Evolved"
NSEOD$Strain<-factor(NSEOD$Strain, levels=c("JCVI-syn1.0","JCVI-syn3B"))


#t-test: Mean of evolved OD versus fixed ancestor value, for syn3B
bNSE<- NSEOD %>% filter(Time=="Evolved") %>% filter(Strain=="JCVI-syn3B")
t.test(x=bNSE$OD,mu=0.079333333,alternative = "g")
#0.001131
#t-test: Mean of evolved OD versus fixed ancestor value, for syn1.0
s1NSE <- NSEOD %>% filter(Time=="Evolved") %>% filter(Strain=="JCVI-syn1.0")
t.test(x=s1NSE$OD,mu=0.876333333,alternative = "g")
#P = 0.0007063

NSEevo<- NSEOD %>% filter(Time=="Evolved")

pgcols<-c("blue","red")
ODfig <- ggplot(NSEevo, aes(x=Strain, y=OD))
ODfig +   geom_errorbar(aes(ymax=c(0.079333333,0.079333333,0.079333333,0.079333333,0.876333333,0.876333333,0.876333333,0.876333333),ymin=c(0.079333333,0.079333333,0.079333333,0.079333333,0.876333333,0.876333333,0.876333333,0.876333333)), color=c("red","red","red","red","blue","blue","blue","blue"),lwd = 2.2,linetype=117) +
  geom_jitter(
  aes(shape = Strain, color = Strain), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = 0.475, dodge.width = 0.2),
  size = 8, stroke=2.3) +
  scale_shape_manual(values = c(0,0)) +
  stat_summary(
    aes(color = Strain),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.5, shape=95,
    position = position_dodge(0.2),
    show.legend=FALSE
  ) +
  scale_color_manual(values = pgcols) +
  labs(x="\nStrain",y="OD 600 after 8 days\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())#+
  #theme(legend.position="none")





mycols<-c("black","black")


pgcols<-c("blue","blue","red","red")
#("red","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue")
ODfig <- ggplot(NSEOD, aes(x=Time, y=OD))
ODfig +  geom_jitter(
    aes(shape = Time, color = Strain),
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.4),
    size = 8, stroke = 2.3, show.legend=FALSE
  ) +
  scale_shape_manual(values = c(0,15)) +
  stat_summary(
    aes(color = c("red","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue")),
    fun.data = "mean_se", fun.args = list(mult = (2)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.8, shape = 95,
    position = position_dodge(0.4),
    show.legend = FALSE) +
  scale_color_manual(values = pgcols) +
  labs(x="\nAssay diet",y= ~ atop(paste("Reproductive fitness (",italic("W'"),")"),###this adds the text that I want
                                  paste(scriptstyle(" ")))) +                           ###however, trying to paste \n messes everything up. Therefore, to gain space between the title and th axis line, I pasted another argument---this starts pn a new line, and, since it was blank, it just functions like \n in thise case
  scale_y_continuous(limits = c(-0.1,2.2), expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())


'''
pgcols<-c("blue","red")
#"red","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"
ODfig <- ggplot(NSEOD, aes(x=Time, y=OD))
ODfig +  geom_jitter(
  aes(shape = Strain, color = Strain),
  position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.4),
  size = 8, stroke = 2.3, show.legend=TRUE
) +
  scale_shape_manual(values = c(15,0)) +
  stat_summary(
    aes(color = Strain),
    fun.data = "mean_se", fun.args = list(mult = (2)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.8, shape = 95,
    position = position_dodge(0.4),
    show.legend = FALSE) +
  scale_color_manual(values = pgcols) +
  labs(x="\nAssay diet",y= ~ atop(paste("Reproductive fitness (",italic("W'"),")"),###this adds the text that I want
                                  paste(scriptstyle(" ")))) +                           ###however, trying to paste \n messes everything up. Therefore, to gain space between the title and th axis line, I pasted another argument---this starts pn a new line, and, since it was blank, it just functions like \n in thise case
  scale_y_continuous(limits = c(-0.1,2.2), expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank()) +
  labs(fill = "TEST")

'''
                                                               
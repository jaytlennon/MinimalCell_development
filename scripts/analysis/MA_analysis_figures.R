rm(list = ls())
#library(ggplot2);library(tidyverse);library(dplyr);library(ggpubr);library(Hmisc)
package.list <- c('vegan', 'ade4', 'viridis', 'gplots', 'BiodiversityR', 'indicspecies', 'ggplot2', 'tidyverse', 'dplyr', 'ggpubr', 'Hmisc')
for (package in package.list) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

myALL<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\MA_summary_for_R.csv")
my3B<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\MA_3B_summary_for_R.csv")
mys1<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\MA_s1_summary_for_R.csv")
MAALL<-as_tibble(myALL)
MA3B<-as_tibble(my3B)
MAs1<-as_tibble(mys1)

MAALL$strain <- as.factor(MAALL$strain)
levels(MAALL$strain)[levels(MAALL$strain)=="MA_s1"] <- "Wildtype"
levels(MAALL$strain)[levels(MAALL$strain)=="MA_3B"] <- "Minimal"
MAALL$strain <- factor(MAALL$strain, levels=c("Wildtype","Minimal"))



####Part 1: test for insertion bias or deletion bias
x1s1<-sum(MAs1$'in')
n1s1<-sum(MAs1$'in')+sum(MAs1$del)
indelprops1 <- prop.test(x = x1s1, n = n1s1, alternative = 'two.sided');indelprops1
#sample prop = 0.215; P-val < 2.2e-16

x13B <- sum(MA3B$'in')
n13B<-sum(MA3B$'in')+sum(MA3B$del)
indelprop3B<-prop.test(x=x13B,n=n13B,p=0.5,alternative="two.sided"); indelprop3B
#sample prop = 0.241; P-val = 4.0e-06

#Significant deletion bias, in terms of number-of-dels > number-of-ins, by chisq test for both strains

twoprop_indel<-prop.test(x=c(x1s1, x13B),n=c(n1s1,n13B),alternative = "two.sided"); twoprop_indel
#p-val > 0.05

#Now test in terms of the total amount of sequence inserted/deleted
x2s1<-sum(MAs1$in_len_tot);n2s1<-x2s1+sum(MAs1$del_len_tot)
indellentotprops1<-prop.test(x=x2s1,n=n2s1,p=0.5,alternative = 'two.sided');indellentotprops1
#sample prop = 0.118, p-val < 2.2e-16
#of the total sequence inserted or deleted, there is a greater proportion deleted than expected by chance, in s1


x23B<-sum(MA3B$in_len_tot);n23B<-x23B+sum(MA3B$del_len_tot)
indellentotprop3B<-prop.test(x=x23B,n=n23B,p=0.5,alternative = 'two.sided');indellentotprop3B
#sample prop = 0.212, p-val <2.2e-16
#As before, conclude that more total seqeunce was deleted, in 3B this time


twoprop_indellentot <- prop.test(x=c(x2s1,x23B),n=c(n2s1,n23B),alternative = 'two.sided');twoprop_indellentot
#p-val < 0.05. Conclude that, in terms of total number of nts inserted or deleted, strain s1 has a greater deletion bias than strain 3B

#Now compare the average sizes of ins among the lines and average sizes of dels within the lines
#test for normality
shapiro.test(na.omit(MAs1$in_len_avg));shapiro.test(na.omit(MAs1$del_len_avg))
ggdensity(na.omit(MAs1$in_len_avg));ggdensity(na.omit(MAs1$del_len_avg))
wilcox.test(x=na.omit(MAs1$in_len_avg),y=na.omit(MAs1$del_len_avg),mu=0,alternative = 'two.sided',paired=F)

#Can I use t-test anyway, since the sample size is large? n>=90. I want to use a t-test because I want to test a hypothesis about a difference in means, not a difference in medians.
#test for equal variances
var.test(x=na.omit(MAs1$in_len_avg),y=na.omit(MAs1$del_len_avg),ratio=1,alternative = 'two.sided')
#F=0.026, p-val < 2.2e-16
indellenavgs1<-t.test(x=na.omit(MAs1$in_len_avg),y=na.omit(MAs1$del_len_avg),mu=0,alternative = 'two.sided',var.equal = F);indellenavgs1
#p-val = 0.0724. For strain s1, marginal deletion bias in terms of avg length of an insertion/deletion.

#test for normality
shapiro.test(na.omit(MA3B$in_len_avg));shapiro.test(na.omit(MA3B$del_len_avg))
ggdensity(na.omit(MA3B$in_len_avg));ggdensity(na.omit(MA3B$del_len_avg))
wilcox.test(x=na.omit(MA3B$in_len_avg),y=na.omit(MA3B$del_len_avg),mu=0,alternative = 'two.sided',paired=F)


#test for equal variances
var.test(x=na.omit(MA3B$in_len_avg),y=na.omit(MA3B$del_len_avg),ratio=1,alternative = 'two.sided')
#F=0.026, p-val < 2.2e-16
indellenavg3B<-t.test(x=na.omit(MA3B$in_len_avg),y=na.omit(MA3B$del_len_avg),mu=0,alternative = 'two.sided',var.equal = F);indellenavg3B
#p-val = 0.333. For strain 3B, no difference in the average size of an in versus a del.




#Test whether the avg len of an insertion differs BETWEEN the strains
shapiro.test(na.omit(MAs1$in_len_avg));shapiro.test(na.omit(MA3B$in_len_avg))
ggdensity(na.omit(MAs1$in_len_avg));ggdensity(na.omit(MA3B$in_len_avg))
wilcox.test(x=na.omit(MAs1$in_len_avg),y=na.omit(MA3B$in_len_avg),mu=0,alternative = 'two.sided',paired=F)

#NO difference in MEDIAN length of an insertion

#test for equal variances
var.test(x=na.omit(MAs1$in_len_avg),y=na.omit(MA3B$in_len_avg),ratio=1,alternative = 'two.sided')
#F=0.208, p-val = 1.132e-05
inlenavgtwostrains<-t.test(x=na.omit(MAs1$in_len_avg),y=na.omit(MA3B$in_len_avg),mu=0,alternative = 'two.sided',var.equal = F);inlenavgtwostrains
#p-val > 0.1. Conclude that the average length of an insertion is not different between the strains


#Test whether the average length of a deletion differs between the strains
shapiro.test(na.omit(MAs1$del_len_avg));shapiro.test(na.omit(MA3B$del_len_avg))
ggdensity(na.omit(MAs1$del_len_avg));ggdensity(na.omit(MA3B$del_len_avg))
wilcox.test(x=na.omit(MAs1$del_len_avg),y=na.omit(MA3B$del_len_avg),mu=0,alternative = 'two.sided',paired=F)
#W=2234.5, p-val = 0.02
#Conclude that the MEDIAN length of a deletion is shorter in 3B

#test for equal variances
var.test(x=na.omit(MAs1$del_len_avg),y=na.omit(MA3B$del_len_avg),ratio=1,alternative = 'two.sided')
#F=38.99, p-val < 2.2 e-16
dellenavgtwostrains<-t.test(x=na.omit(MAs1$del_len_avg),y=na.omit(MA3B$del_len_avg),mu=0,alternative = 'two.sided',var.equal = F);dellenavgtwostrains
#p-val = 0.130. Conclude that the average length of a deletion is not different between the strains




##########Length of a deletion:
#What if we include the ones longer than 100 bp?
var.test(x=na.omit(MAs1$ALLdel_len_avg),y=na.omit(MA3B$ALLdel_len_avg),ratio=1,alternative = 'two.sided')
#F=8486, p-val < 2.2 e-16
ALLdellenavgtwostrains<-t.test(x=na.omit(MAs1$ALLdel_len_avg),y=na.omit(MA3B$ALLdel_len_avg),mu=0,alternative = 'two.sided',var.equal = F);ALLdellenavgtwostrains
#t = 2.18, p=0.048. Conclude that if you include ALL deletions, the mean deletion length is longer in syn1.0
wilcox.test(x=na.omit(MAs1$ALLdel_len_avg),y=na.omit(MA3B$ALLdel_len_avg),mu=0,alternative = 'two.sided',paired=F)
#The median length is of a deletion is marginally longer in s1, with p < 0.1

####
#Test whether there the avg length of a deletion is greater than the length of an insertion, including ALL deletions, for strain s1
shapiro.test(na.omit(MAs1$ALLin_len_avg));shapiro.test(na.omit(MAs1$ALLdel_len_avg))
ggdensity(na.omit(MAs1$ALLin_len_avg));ggdensity(na.omit(MAs1$ALLdel_len_avg))
wilcox.test(x=na.omit(MAs1$ALLin_len_avg),y=na.omit(MAs1$ALLdel_len_avg),mu=0,alternative = 'two.sided',paired=F)

#Can I use t-test anyway, since the sample size is large? n>=90. I want to use a t-test because I want to test a hypothesis about a difference in means, not a difference in medians.
#test for equal variances
var.test(x=na.omit(MAs1$ALLin_len_avg),y=na.omit(MAs1$ALLdel_len_avg),ratio=1,alternative = 'two.sided')
#F=0.026, p-val < 2.2e-16
ALLindellenavgs1<-t.test(x=na.omit(MAs1$ALLin_len_avg),y=na.omit(MAs1$ALLdel_len_avg),mu=0,alternative = 'two.sided',var.equal = F);ALLindellenavgs1
#p-val = 0.0459. Now when ALL are included, for strain s1, SIGNIFICANT deletion bias in terms of avg length of an insertion/deletion.




#Is the avg length of a deletion longer in the wildtype than in the minimal?
shapiro.test(na.omit(MAs1$ALLdel_len_avg));shapiro.test(na.omit(MA3B$ALLdel_len_avg))
ggdensity(na.omit(MAs1$ALLdel_len_avg));ggdensity(na.omit(MA3B$ALLdel_len_avg))
wilcox.test(x=na.omit(MAs1$ALLdel_len_avg),y=na.omit(MA3B$ALLdel_len_avg),mu=0,alternative = 'two.sided',paired=F)
############
var.test(x=na.omit(MAs1$ALLdel_len_avg),y=na.omit(MA3B$ALLdel_len_avg),ratio=1,alternative = 'two.sided')
#F=0.026, p-val < 2.2e-16
ALLindellenavgs1<-t.test(x=na.omit(MAs1$ALLdel_len_avg),y=na.omit(MA3B$ALLdel_len_avg),mu=0,alternative = 'two.sided',var.equal = F);ALLindellenavgs1


#repeat the indel t-test analyses in a 2-way anova framework, using ALL deletions
#myindellenavgs<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\MA_indel_avg_len.csv")
myindellenavgs<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\MA_ALLindel_avg_len_tools.csv")
MAindels<-as_tibble(na.omit(myindellenavgs))
indel2way<-aov(MAindels$len_avg ~ MAindels$strain + MAindels$type + MAindels$strain*MAindels$type)
summary(indel2way)
TukeyHSD(indel2way)
#no significant differences is this naive framework


###############################
#Brief detour: Repeat the histograms of insertion length/deletion length, but including indels larger than 100 bp
#shapiro.test(na.omit(MAs1$ALLin_len_avg));#shapiro.test(na.omit(MA3b$indel_len_avg))
ggdensity(na.omit(MAs1$ALLin_len_avg))
ggdensity(na.omit(MA3B$ALLin_len_avg))
ggdensity(na.omit(MAs1$ALLdel_len_avg))
ggdensity(na.omit(MA3B$ALLdel_len_avg))


####PART 2: investigate to-AT-biased or to-GC-biased mutation spectrum
x3s1<-sum(MAs1$to_AT_tot);n3s1<-x3s1+sum(MAs1$to_CG_tot)
#toATprops1<-prop.test(x=x3s1,n=n3s1,p=(1-0.760011272),alternative = 'two.sided');toATprops1
toATprops1<-prop.test(x=x3s1,n=n3s1,p=(1-0.76),alternative = 'two.sided');toATprops1
#sample prop = 0.903, p-val < 2.2e-16

x33B<-sum(MA3B$to_AT_tot);n33B<-x33B+sum(MA3B$to_CG_tot)
#toATprop3B<-prop.test(x=x33B,n=n33B,p=(1-0.757152926),alternative = 'two.sided');toATprop3B
toATprop3B<-prop.test(x=x33B,n=n33B,p=(1-0.76),alternative = 'two.sided');toATprop3B
#sample prop = 0.970, p-val <2.2e-16

#the two strains have almost identical expectation proportions. I will go ahead and compare them.
twoprop_ATbias<-prop.test(x=c(x3s1, x33B),n=c(n3s1,n33B),alternative = "two.sided"); twoprop_ATbias
#p-val = 3.1e-06 signif difference suggesting that AT bias is stronger in syn3B. But also syn3B had a higher expectation proportion which I wasnt able to account for----not sure what to do. We do know that they are signifly different because the confidence interval does not include the actual expected difference--- CI = [-0.089, -0.045], when the expected difference is -0.003
#Expected difference is -0.002858346 (i.e., higher for 3B).
#95 percent confidence interval:
# -0.08939673 -0.04458173
#The confidence interval does not include the expected difference, which is -0.003.
#conclude: Both strains have AT biased mutation. The AT bias is stronger in 3B.

#there was an excess of AT mutation bias in the minimal cell. Was this reflected in the mutation spectrum of mutations in significant genes in the NSE?
x3s1_NSE<-23;n3s1_NSE<-23+9
toATprops1_NSE<-prop.test(x=x3s1_NSE,n=n3s1_NSE,p=(1-0.76),alternative = 'two.sided');toATprops1_NSE

x33B_NSE<-34;n33B_NSE<-34+1
#toATprop3B<-prop.test(x=x33B,n=n33B,p=(1-0.757152926),alternative = 'two.sided');toATprop3B
toATprop3B_NSE<-prop.test(x=x33B_NSE,n=n33B_NSE,p=(1-0.76),alternative = 'two.sided');toATprop3B_NSE
#stronger AT bias in the NSE's significant genes

twoprop_ATbias_NSE<-prop.test(x=c(x3s1_NSE, x33B_NSE),n=c(n3s1_NSE,n33B_NSE),alternative = "two.sided"); twoprop_ATbias_NSE

##the gene ung was deleted from the minimal cell. Is there a higher proportion of C-->T mutations in the minimal cell's MA? In its significant NSE genes?
x3s1_MA_ung<-759;n3s1_MA_ung<-1591
x33B_MA_ung<-341;n33B_MA_ung<-521
twoprop_ATbias_MA_ung<-prop.test(x=c(x3s1_MA_ung, x33B_MA_ung),n=c(n3s1_MA_ung,n33B_MA_ung),alternative = "two.sided"); twoprop_ATbias_MA_ung

x3s1_NSE_ung<-18;n3s1_NSE_ung<-41
x33B_NSE_ung<-21;n33B_NSE_ung<-43
twoprop_ATbias_NSE_ung<-prop.test(x=c(x3s1_NSE_ung, x33B_NSE_ung),n=c(n3s1_NSE_ung,n33B_NSE_ung),alternative = "two.sided"); twoprop_ATbias_NSE_ung
#no difference in the relative proportions of C-->T mutations between the strains

#PART 3: transition/transversion ratios
x4s1<-sum(MAs1$ts);n4s1<-x4s1+sum(MAs1$tv)
tstvs1<-prop.test(x=x4s1,n=n4s1,p=(1/3),alternative = 'two.sided');tstvs1
#sample prop = 0.492, p-val < 2.2e-16
#conclude significant transition bias in s1

x43B<-sum(MA3B$ts);n43B<-x43B+sum(MA3B$tv)
tstv3B<-prop.test(x=x43B,n=n43B,p=(1/3),alternative = 'two.sided');tstv3B
#sample prop = 0.670, p-val <2.2e-16
#conclude signif transition bias in 3B

twoprop_tstv <- prop.test(x=c(x4s1,x43B),n=c(n4s1,n43B),alternative = 'two.sided');twoprop_tstv
#p-val < 1.0e-12. Conclude that transition bias is significantly stronger in 3B




#PART 4: Do a greater proportion of lines possess a large (> 100 bp) mutation in one strain versus the other?
max(MAs1$SV);max(MA3B$SV)
#since the max is 1, I can just take a sum to get the number of lines with a large mut
x5s1=sum(MAs1$SV);n5s1=90;x53B=sum(MA3B$SV);n53B=67
twoprop_large <- prop.test(x=c(x5s1,x53B),n=c(n5s1,n43B),alternative = 'two.sided');twoprop_large
#prop 1: 0.188888888888; prop 2: 0; p-val < e-09




#PART 5: Compare whether the total mutn rate was different between the strains

#first, per nt (which is the interesting comparison as far as cell biology goes)
#test for normality
shapiro.test(MAs1$per_nt_per_gen);shapiro.test(MA3B$per_nt_per_gen)
ggdensity(MAs1$per_nt_per_gen);ggdensity(MA3B$per_nt_per_gen)
wilcox.test(x=MAs1$per_nt_per_gen,y=MA3B$per_nt_per_gen,mu=0,alternative = 'two.sided',paired=F)
#P > 0.05

shapiro.test(log(MAs1$per_nt_per_gen));shapiro.test(log(MA3B$per_nt_per_gen))
ggdensity(log(MAs1$per_nt_per_gen));ggdensity(log(MA3B$per_nt_per_gen))
wilcox.test(x=log(MAs1$per_nt_per_gen),y=log(MA3B$per_nt_per_gen),mu=0,alternative = 'two.sided',paired=F)
#W=same--- log transform doesnt affect the median


#Can I use t-test anyway, since the sample size is large? n>=90. I want to use a t-test because I want to test a hypothesis about a difference in means, not a difference in medians.
#test for equal variances
var.test(x=MAs1$per_nt_per_gen,y=MA3B$per_nt_per_gen,ratio=1,alternative = 'two.sided')
#p >0.05
pernt_twostrains<-t.test(x=MAs1$per_nt_per_gen,y=MA3B$per_nt_per_gen,mu=0,alternative = 'two.sided',var.equal = T);pernt_twostrains
#p > 0.05, = 0.541

var.test(x=log(MAs1$per_nt_per_gen),y=log(MA3B$per_nt_per_gen),ratio=1,alternative = 'two.sided')
#F=0.717, p-val = 0.144
pernt_twostrains_log<-t.test(x=log(MAs1$per_nt_per_gen),y=log(MA3B$per_nt_per_gen),mu=0,alternative = 'two.sided',var.equal = T);pernt_twostrains_log
#p-val = .65
#Conclude that mutation rate is not significantly higher in the 3B strain. So this suggests there was no great effect of removing those extra DNA replication/repair genes. 

#######
###Before moving on: Let's compare the bps-only mutation rates.
shapiro.test(log(MAs1$bpsub_per_nt));shapiro.test(log(MA3B$bpsub_per_nt))
ggdensity(log(MAs1$bpsub_per_nt));ggdensity(log(MA3B$bpsub_per_nt))
wilcox.test(x=log(MAs1$bpsub_per_nt),y=log(MA3B$bpsub_per_nt),mu=0,alternative = 'two.sided',paired=F)
#W=2409,p=0.957
var.test(x=log(MAs1$bpsub_per_nt),y=log(MA3B$bpsub_per_nt),ratio=1,alternative = 'two.sided')
#F=0.852, p-val = 0.500
bpsub_pernt_twostrains_log<-t.test(x=log(MAs1$bpsub_per_nt),y=log(MA3B$bpsub_per_nt),mu=0,alternative = 'two.sided',var.equal = T);bpsub_pernt_twostrains_log


#t=0.099, p-val = 0.921
###Conclude that there is no different in the bps rate
var.test(x=(MAs1$bpsub_per_nt),y=(MA3B$bpsub_per_nt),ratio=1,alternative = 'two.sided')
bpsub_pernt_twostrains<-t.test(x=(MAs1$bpsub_per_nt),y=(MA3B$bpsub_per_nt),mu=0,alternative = 'two.sided',var.equal = T);bpsub_pernt_twostrains
#t=0.051, p-val = 0.959. Same conclusion. No difference in the bps rate.


#What COULD have affecte the NSE, is that syn3B prolly had a lower per genome mutation rate---that's the actual input of mutations that could drive evolution.


#second, per total CDS. This is Mike Lynch's evolutionary question---a test of the DBH---and DBH is not testable here, because the minimal cell syn3B has not had any evolution.
#HOWEVER it is an interesting question to the extent that it would have altered the dynamics during the NSE!
#test for normality
shapiro.test(MAs1$per_CDS_per_gen);shapiro.test(MA3B$per_CDS_per_gen)
ggdensity(MAs1$per_CDS_per_gen);ggdensity(MA3B$per_CDS_per_gen)
wilcox.test(x=MAs1$per_CDS_per_gen,y=MA3B$per_CDS_per_gen,mu=0,alternative = 'two.sided',paired=F)
#W=6450, p-val = 6.484e-09

shapiro.test(log(MAs1$per_CDS_per_gen));shapiro.test(log(MA3B$per_CDS_per_gen))
ggdensity(log(MAs1$per_CDS_per_gen));ggdensity(log(MA3B$per_CDS_per_gen))
wilcox.test(x=log(MAs1$per_CDS_per_gen),y=log(MA3B$per_CDS_per_gen),mu=0,alternative = 'two.sided',paired=F)
#W=6450, p-val = 6.484e-09

#test for equal variances
var.test(x=MAs1$per_CDS_per_gen,y=MA3B$per_CDS_per_gen,ratio=1,alternative = 'two.sided')
#F=1.560, p-val = 0.0337
perCDS_twostrains<-t.test(x=MAs1$per_CDS_per_gen,y=MA3B$per_CDS_per_gen,mu=0,alternative = 'two.sided',var.equal = F);perCDS_twostrains
#p-val = 8.699e-09

var.test(x=log(MAs1$per_CDS_per_gen),y=log(MA3B$per_CDS_per_gen),ratio=1,alternative = 'two.sided')
#F=0.816, p-val = 0.334
perCDS_twostrains_log<-t.test(x=log(MAs1$per_CDS_per_gen),y=log(MA3B$per_CDS_per_gen),mu=0,alternative = 'two.sided',var.equal = T);perCDS_twostrains_log
#p-val < 0.05




#ok so part 2b is per GENOMEper gen. This is more impt than per CDS for the NSE I think (even tho per CDS is more important for Mike Lynch's comparative question!)
shapiro.test(MAs1$per_genome_per_gen);shapiro.test(MA3B$per_genome_per_gen)
ggdensity(MAs1$per_genome_per_gen);ggdensity(MA3B$per_genome_per_gen)
wilcox.test(x=MAs1$per_genome_per_gen,y=MA3B$per_genome_per_gen,mu=0,alternative = 'two.sided',paired=F)
#W=6450, p-val = 6.484e-09

shapiro.test(log(MAs1$per_genome_per_gen));shapiro.test(log(MA3B$per_genome_per_gen))
ggdensity(log(MAs1$per_genome_per_gen));ggdensity(log(MA3B$per_genome_per_gen))
wilcox.test(x=log(MAs1$per_genome_per_gen),y=log(MA3B$per_genome_per_gen),mu=0,alternative = 'two.sided',paired=F)
#W=6450, p-val = 6.484e-09

#test for equal variances
var.test(x=MAs1$per_genome_per_gen,y=MA3B$per_genome_per_gen,ratio=1,alternative = 'two.sided')
#F=1.649, p-val = 0.0170
pergenome_twostrains<-t.test(x=MAs1$per_genome_per_gen,y=MA3B$per_genome_per_gen,mu=0,alternative = 'two.sided',var.equal = F);pergenome_twostrains
#p-val = 4.924e-10

var.test(x=log(MAs1$per_genome_per_gen),y=log(MA3B$per_genome_per_gen),ratio=1,alternative = 'two.sided')
#F=0.816, p-val = 0.334
pergenome_twostrains_log<-t.test(x=log(MAs1$per_genome_per_gen),y=log(MA3B$per_genome_per_gen),mu=0,alternative = 'two.sided',var.equal = T);pergenome_twostrains_log
#p-val < 2.2e-16


#conclude that the minimal cell has a significantly lower per-genome-per-gen mutation rate. Is a difference of 0.0338 versus 0.0177 mutations per genome per generation biologically signif? In the context of the NSE, this would result in a difference of 
0.033790614*2000*30000000
0.017655872*2000*30000000
#2.0e09 vs 1.1e09 mutations tested during the course of the evolution, assuming an avg cell dens of 1e07 cells/mL. So either way, a mutation at every position has been tested 2.0e09/1e06 = 2000 (wildtype) versus 1.1e09/5e05=2200 (minimal) times during the NSE.

############################################################
#Part 6: i) Test dN/dS ratio. ii) Compare proportion of called SNMs in CDS regions versus non-CDS regions to the theoretical null expectation.

#First look at dN/dS value for each line
shapiro.test(MAs1$dN_dS);shapiro.test(MA3B$dN_dS)
#non-normal distrns
wilcox.test(x=MAs1$dN_dS,y=NULL,mu=1,alternative="two.sided")
#P=1.703e-10---- dN/dS is significantly different (greater) than 1
wilcox.test(x=MA3B$dN_dS,y=NULL,mu=1,alternative="two.sided")
#P=0.02665----- dN/dS is significantly different (greater) than 1
t.test(x=MAs1$dN_dS,y=NULL,mu=1,alternative="two.sided")
t.test(x=MA3B$dN_dS,y=NULL,mu=1,alternative="two.sided")
#Same inference when using t-test

#Next, for each strain, look at the total number of CDS SNMs and non-CDS SNMs. Do they occur in proportions significantly different from the null expectation?
x6s1<-sum(MAs1$in_CDS)
n6s1<-sum(MAs1$in_CDS)+sum(MAs1$not_in_CDS)
CDSprops1 <- prop.test(x = x6s1, n = n6s1, p = (926310/1078809), alternative = 'two.sided');CDSprops1
#P = 0.876. There is not significant deviation from the proportion of mutations falling in CDS regions expected due to chance alone. 86% compared to expected 86%

x63B<-sum(MA3B$in_CDS)
n63B<-sum(MA3B$in_CDS)+sum(MA3B$not_in_CDS)
CDSprops1 <- prop.test(x = x63B, n = n63B, p = (479721/543379), alternative = 'two.sided');CDSprops1
#P 0.069. There is not significant deviation from the proportion of mutations falling in CDS regions expected due to chance alone. 86%, compared to expected 88%



######################
#Make some figures

mySNM<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\MA_grouped_bar_SNM.csv")
mymutp<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\MA_grouped_bar_mut.types.csv")
MASNM<-as_tibble(mySNM)
MAmutp<-as_tibble(mymutp)
NSESNM<-as_tibble(myNSE_SNM)

MASNM$strain <- as.factor(MASNM$strain)
MASNM$mut <- as.factor(MASNM$mut)
levels(MASNM$strain)[levels(MASNM$strain)=="s1"] <- "Wildtype"
levels(MASNM$strain)[levels(MASNM$strain)=="3B"] <- "Minimal"
MASNM$strain <- factor(MASNM$strain, levels=c("Wildtype","Minimal"))
MASNM$mut <- factor(MASNM$mut, levels = c("A:T to C:G","A:T to G:C","A:T to T:A","C:G to G:C", "C:G to T:A", "C:G to A:T"))

myNSE_SNM<-read_csv("C:\\Users\\rmoge\\Box Sync\\Mycoplasma\\Experiments\\data\\MA\\NSE_grouped_bar_SNM.csv")
NSESNM$strain <- as.factor(NSESNM$strain)
NSESNM$mut <- as.factor(NSESNM$mut)
levels(NSESNM$strain)[levels(NSESNM$strain)=="s1"] <- "Wildtype"
levels(NSESNM$strain)[levels(NSESNM$strain)=="3B"] <- "Minimal"
NSESNM$strain <- factor(NSESNM$strain, levels=c("Wildtype","Minimal"))
NSESNM$mut <- factor(NSESNM$mut, levels = c("A:T to C:G","A:T to G:C","A:T to T:A","C:G to G:C", "C:G to T:A", "C:G to A:T"))

MAmutp$strain <- as.factor(MAmutp$strain)
MAmutp$mut <- as.factor(MAmutp$mut)
levels(MAmutp$strain)[levels(MAmutp$strain)=="s1"] <- "Wildtype"
levels(MAmutp$strain)[levels(MAmutp$strain)=="3B"] <- "Minimal"
levels(MAmutp$mut)[levels(MAmutp$mut)=="SNM"] <- "Single-nucleotide"
MAmutp$strain <- factor(MAmutp$strain, levels=c("Wildtype","Minimal"))
MAmutp$mut <- factor(MAmutp$mut, levels = c("Insertion","Deletion","Over 100 bp","Single-nucleotide"))

MASNM<-MASNM %>% mutate(prop2 = sprintf("%0.2f",prop))
NSESNM<-NSESNM %>% mutate(prop2 = sprintf("%0.2f",prop))
MAmutp<-MAmutp %>% mutate(prop2 = sprintf("%0.2f",prop))

#theme_set(theme_pubclean())


ggplot(MASNM, aes(x = mut, y = prop)) +
  geom_bar(
    aes(color = strain, fill = strain),
    stat = "identity", position = position_stack()
  ) +
  scale_color_manual(values = c("blue", "red"))+
  scale_fill_manual(values = c("blue", "red"))

# Use position = position_dodge() 
p <- ggplot(MASNM, aes(x = mut, y = prop)) +
  geom_bar(
    aes(color = strain, fill = strain),
    stat = "identity", position = position_dodge(0.5),
    width = 0.35
  ) +
  scale_color_manual(values = c("blue", "red"))+
  scale_fill_manual(values = c("blue", "red")) +
  geom_text(aes(label = prop2, group = strain), position = position_dodge(0.5), vjust = -0.3, size = 5.5)+
  labs(x="\nMutation type",y="Proportion\n") +
  theme(axis.text.x = element_text(angle=-45, vjust=0,hjust=0.2, size =26))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit(.3, "cm"), axis.ticks.x = element_blank())
p

ggplot(NSESNM, aes(x = mut, y = prop)) +
  geom_bar(
    aes(color = strain, fill = strain),
    stat = "identity", position = position_stack()
  ) +
  scale_color_manual(values = c("blue", "red"))+
  scale_fill_manual(values = c("blue", "red"))
# Use position = position_dodge() 
NSEp <- ggplot(NSESNM, aes(x = mut, y = prop)) +
  geom_bar(
    aes(color = strain, fill = strain),
    stat = "identity", position = position_dodge(0.5),
    width = 0.3
  ) +
  scale_color_manual(values = c("blue", "red"))+
  scale_fill_manual(values = c("blue", "red")) +
  geom_text(aes(label = prop2, group = strain), position = position_dodge(0.5), vjust = -0.3, size = 5.5)+
  labs(x="\nMutation type",y="Proportion\n") +
  theme(axis.text.x = element_text(angle=-45, vjust=0,hjust=0.2, size =26))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit(.3, "cm"), axis.ticks.x = element_blank())
NSEp


ggplot(MAmutp, aes(x = mut, y = prop)) +
  geom_bar(
    aes(color = strain, fill = strain),
    stat = "identity", position = position_stack()
  ) +
  scale_color_manual(values = c("blue", "red"))+
  scale_fill_manual(values = c("blue", "red"))

# Use position = position_dodge() 
p <- ggplot(MAmutp, aes(x = mut, y = prop)) +
  geom_bar(
    aes(color = strain, fill = strain),
    stat = "identity", position = position_dodge(0.5),
    width = 0.35
  ) +
  scale_color_manual(values = c("blue", "red"))+
  scale_fill_manual(values = c("blue", "red")) +
  geom_text(aes(label = prop2, group = strain), position = position_dodge(0.5), vjust = -0.3, size = 5.5)+
  labs(x="\nMutation type",y="Proportion\n") +
  theme(axis.text.x = element_text(angle=-45, vjust=0,hjust=0.2, size =26))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),legend.key.size = unit(2, "lines"),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit(.3, "cm"), axis.ticks.x = element_blank())
p


#Now do two strip charts, comparing the mutation rate, view thru the per nt lens and thru the per CDS lens
mycols<-c('blue','red')
pernt <- ggplot(MAALL, aes(x=strain, y=per_nt_per_gen)) +
  geom_jitter(
  aes(shape = strain, color = strain), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
  position = position_jitterdodge(jitter.width = .9, dodge.width = .2),
  size = 1, stroke = 2.3) +
  scale_shape_manual(values = c(15,15)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = strain),
    fun.data = "mean_se", fun.args = list(mult = (2)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.5, shape=95,
    position = position_dodge(0.2),
    show.legend = FALSE
  ) +
  scale_color_manual(values = mycols) +
  labs(x="\nStrain",y="muts/nt/gen\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",legend.key=element_blank(),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line = element_line(colour = "black"), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())
pernt


#Now do two strip charts, comparing the mutation rate, view thru the per nt lens and thru the per CDS lens
mycols<-c('blue','red')
perCDS <- ggplot(MAALL, aes(x=strain, y=per_CDS_per_gen)) +
  geom_jitter(
    aes(shape = strain, color = strain), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
    position = position_jitterdodge(jitter.width = .9, dodge.width = .2),
    size = 1, stroke = 2.3) +
  scale_shape_manual(values = c(15,15)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = strain),
    fun.data = "mean_se", fun.args = list(mult = (2)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.5, shape=95,
    position = position_dodge(0.2),
    show.legend = FALSE
  ) +
  scale_color_manual(values = mycols) +
  labs(x="\nStrain",y="muts/CDS/gen\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",legend.key=element_blank(),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line = element_line(colour = "black"), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())
perCDS

mycols<-c('blue','red')
pergenome <- ggplot(MAALL, aes(x=strain, y=per_genome_per_gen)) +
  geom_jitter(
    aes(shape = strain, color = strain), #I cut out the "shape = Evolution, " part of the aes so that all shapes would be the same shape
    position = position_jitterdodge(jitter.width = .9, dodge.width = .2),
    size = 1, stroke = 2.3) +
  scale_shape_manual(values = c(15,15)) +#12 is a square with a vertical cross inside it
  stat_summary(
    aes(color = strain),
    fun.data = "mean_se", fun.args = list(mult = (2)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 1.5, shape=95,
    position = position_dodge(0.2),
    show.legend = FALSE
  ) +
  scale_color_manual(values = mycols) +
  labs(x="\nStrain",y="muts/genome/gen\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",legend.key=element_blank(),axis.text=element_text(size=34),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.line = element_line(colour = "black"), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_line(color = "black", size = 3.5), axis.ticks.length = unit (.3, "cm"), axis.ticks.x = element_blank())
pergenome

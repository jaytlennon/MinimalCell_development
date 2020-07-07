# -*- coding: utf-8 -*-
"""
Created on 2020/02/19

@author: Roy Moger-Reischer
"""

#%%
from __future__ import division
import re, os,sys, math, operator,random, copy,collections,time; import numpy as np; import pandas as pd; import csv
from itertools import groupby; import pprint as pp
import matplotlib.pyplot as plt
#%%
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
#%%
import vcf
#%%
#mytest=pd.read_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_testing\3B_1_clean.VCF')

mytest = vcf.Reader(open(r'C:\Users\rmoge\GitHub\MuriSam\MA_testing\3B_1_clean.VCF','r'))
for r in mytest:
    print(r)
    
    
#%%
#Let's make a matrix of ALL LINES by sites mutated---it's a matter of 0s and 1s
def get_big_mut_matrix(infile="",prefix="3B_"):
    straindict=dict()
    for n in range(1,97):
        try:
            straindict[n]=dict()
            stringname=prefix + str(n)
            thisfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_' + prefix + 'VCF\\' + prefix + str(n) + r'_clean.VCF'
            myvcf = vcf.Reader(open(thisfile,'r'))
            for m in myvcf:
                straindict[n][m.POS]=1
        except FileNotFoundError:
            continue
                
            
    return straindict


mutlist_3B=get_big_mut_matrix(prefix="3B_")
MA_3B_mutpd=pd.DataFrame.from_dict(mutlist_3B,orient="index")
MA_3B_mutpd.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_3B_annotatedGD\MA_3B_all.muts.csv')
        
mutlist_s1=get_big_mut_matrix(prefix="s1_")
MA_s1_mutpd=pd.DataFrame.from_dict(mutlist_s1,orient="index")
MA_s1_mutpd.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_annotatedGD\MA_s1_all.muts.csv')        

#%%
def count_fixed_mutations_3B(infile="",samplename=""):
    sample = vcf.Reader(open(infile,'r'))
    basis=dict()
    basis['in']=0;basis['del']=0;basis["SV"]=0;basis["AT_CG"]=0;basis["AT_GC"]=0;basis["AT_TA"]=0;basis["CG_GC"]=0;basis["CG_TA"]=0;basis["GC_TA"]=0
    basis["SVdel"]=0;basis["SVin"]=0
    del_lens=list();in_lens=list();SVdel_lens=list();SVin_lens=list();bothdelslens=list();bothinslens=list()
    for m in sample:
        if m.POS >= 76900 and m.POS <= 76909:#there is an oligo-T deletion in the MA_3B ancestor, but breseq doesn't always score it exactly the same, and therefore it wasn't subtracted out of the lines---hence I am subtracting it here
            continue
        elif m.POS == 20938  or m.POS==484801 or m.POS==484793 or m.POS==484795 or m.POS==516354 or m.POS==484791 or m.POS==516356 or m.POS==516362 or m.POS==420163 or m.POS==20023 or m.POS==14906 or m.POS==385785 or m.POS==76905 or m.POS==145853 or m.POS==437910 or m.POS==437909 or m.POS==449149:
            continue
        else:
            if m.INFO["AF"][0]!=1.0:
                #print("not fixed")
                continue
            else:
                #print(m.ALT[0])
                if len(m.ALT[0])!=1 or len(m.REF)!=1:
                    #figure out if it's an in, a del, or a SV
                    if len(m.ALT[0])>100 or len(m.REF)>100:
                        basis["SV"]+=1
                        if len(m.ALT[0]) > len(m.REF):
                            basis["SVin"]+=1
                            SVin_lens.append(len(m.ALT[0]) - len(m.REF))
                            bothinslens.append(len(m.ALT[0]) - len(m.REF))
                        elif len(m.ALT[0]) < len(m.REF):
                            basis["SVdel"]+=1
                            SVdel_lens.append(len(m.REF) - len(m.ALT[0]))
                            bothdelslens.append(len(m.REF) - len(m.ALT[0]))
                    elif len(m.ALT[0]) > len(m.REF):
                        basis['in']+=1
                        in_lens.append(len(m.ALT[0]) - len(m.REF))
                        bothinslens.append(len(m.ALT[0]) - len(m.REF))
                    elif len(m.ALT[0]) < len(m.REF):
                        basis['del']+=1
                        del_lens.append(len(m.REF) - len(m.ALT[0]))
                        bothdelslens.append(len(m.REF) - len(m.ALT[0]))
                    else:
                        continue
                else:
                    if m.REF == "A":
                        if m.ALT[0] == "C":
                            basis["AT_CG"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_GC"]+=1
                        else:
                            basis["AT_TA"]+=1
                    elif m.REF == "C":
                        if m.ALT[0] == "A":
                            basis["GC_TA"]+=1
                        elif m.ALT[0] == "G":
                            basis["CG_GC"]+=1
                        else:
                            basis["CG_TA"]+=1
                    elif m.REF == "G":
                        if m.ALT[0] == "A":
                            basis["CG_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["CG_GC"]+=1
                        elif m.ALT[0] == "T":
                            basis["GC_TA"]+=1
                    else:
                        if m.ALT[0] == "A":
                            basis["AT_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["AT_GC"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_CG"]
                    #figure out which of the 6 possible substs it was
    '''if len(in_lens)==0:
        in_lens.append(0)
    else:
        pass
    if len(del_lens)==0:
        del_lens.append(0)
    else:
        pass'''
    basis['in_len_avg']=np.average(in_lens);basis['del_len_avg']=np.average(del_lens)
    basis['in_len_tot']=np.sum(in_lens);basis['del_len_tot']=np.sum(del_lens)
    basis['SVin_len_avg']=np.average(SVin_lens);basis['SVdel_len_avg']=np.average(SVdel_lens)
    basis['SVin_len_tot']=np.sum(in_lens);basis['SVdel_len_tot']=np.sum(del_lens)
    basis['ALLin_len_avg']=np.average(bothinslens);basis['ALLdel_len_avg']=np.average(bothdelslens)
    basis['ALLin_len_tot']=np.sum(bothinslens);basis['ALLdel_len_tot']=np.sum(bothdelslens)
    return basis
#%%
count_fixed_mutations_3B(r'C:\Users\rmoge\GitHub\MuriSam\MA_testing\3B_1_clean.VCF')

#%%
wholedict=dict()
for n in range(1,97):
    stringname="3B_" + str(n)
    try:
        thisfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_3B_VCF\3B_' + str(n) + r'_clean.VCF'
        wholedict[stringname]=count_fixed_mutations_3B(thisfile)
    except FileNotFoundError:
        continue        
MA_3B_df=pd.DataFrame.from_dict(wholedict,orient="index")
MA_3B_df.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_3B_VCF\MA_3B_summary.csv')


#%%
def count_fixed_mutations_s1(infile="",samplename=""):
    sample = vcf.Reader(open(infile,'r'))
    basis=dict()
    basis['in']=0;basis['del']=0;basis["SV"]=0;basis["AT_CG"]=0;basis["AT_GC"]=0;basis["AT_TA"]=0;basis["CG_GC"]=0;basis["CG_TA"]=0;basis["GC_TA"]=0
    basis["SVdel"]=0;basis["SVin"]=0
    del_lens=list();in_lens=list();SVdel_lens=list();SVin_lens=list();bothdelslens=list();bothinslens=list()
    for m in sample:
        if m.POS == 917603 or m.POS==120719 or m.POS==359605 or m.POS==549055 or m.POS==118607 or m.POS==538105 or m.POS==855298:#there is mutation at this position in the ancestor that was not caught by gdtools subtract because breseq didnt catch it
            continue
        elif m.POS == 31655:
            continue
        else:
            if m.INFO["AF"][0]!=1.0:
                #print("not fixed")
                continue
            else:
                #print(m.ALT[0])
                if len(m.ALT[0])!=1 or len(m.REF)!=1:
                    #figure out if it's an in, a del, or a SV
                    if len(m.ALT[0])>100 or len(m.REF)>100:
                        basis["SV"]+=1
                        if len(m.ALT[0]) > len(m.REF):
                            basis["SVin"]+=1
                            SVin_lens.append(len(m.ALT[0]) - len(m.REF))
                            bothinslens.append(len(m.ALT[0]) - len(m.REF))
                        elif len(m.ALT[0]) < len(m.REF):
                            basis["SVdel"]+=1
                            SVdel_lens.append(len(m.REF) - len(m.ALT[0]))
                            bothdelslens.append(len(m.REF) - len(m.ALT[0]))
                    elif len(m.ALT[0]) > len(m.REF):
                        basis['in']+=1
                        in_lens.append(len(m.ALT[0]) - len(m.REF))
                        bothinslens.append(len(m.ALT[0]) - len(m.REF))
                    elif len(m.ALT[0]) < len(m.REF):
                        basis['del']+=1
                        del_lens.append(len(m.REF) - len(m.ALT[0]))
                        bothdelslens.append(len(m.REF) - len(m.ALT[0]))
                    else:
                        continue
                else:
                    if m.REF == "A":
                        if m.ALT[0] == "C":
                            basis["AT_CG"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_GC"]+=1
                        else:
                            basis["AT_TA"]+=1
                    elif m.REF == "C":
                        if m.ALT[0] == "A":
                            basis["GC_TA"]+=1
                        elif m.ALT[0] == "G":
                            basis["CG_GC"]+=1
                        else:
                            basis["CG_TA"]+=1
                    elif m.REF == "G":
                        if m.ALT[0] == "A":
                            basis["CG_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["CG_GC"]+=1
                        elif m.ALT[0] == "T":
                            basis["GC_TA"]+=1
                    else:
                        if m.ALT[0] == "A":
                            basis["AT_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["AT_GC"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_CG"]
                    #figure out which of the 6 possible substs it was
    '''if len(in_lens)==0:
        in_lens.append(0)
    else:
        pass
    if len(del_lens)==0:
        del_lens.append(0)
    else:
        pass'''
    basis['in_len_avg']=np.average(in_lens);basis['del_len_avg']=np.average(del_lens)
    basis['in_len_tot']=np.sum(in_lens);basis['del_len_tot']=np.sum(del_lens)
    basis['SVin_len_avg']=np.average(SVin_lens);basis['SVdel_len_avg']=np.average(SVdel_lens)
    basis['SVin_len_tot']=np.sum(in_lens);basis['SVdel_len_tot']=np.sum(del_lens)
    basis['ALLin_len_avg']=np.average(bothinslens);basis['ALLdel_len_avg']=np.average(bothdelslens)
    basis['ALLin_len_tot']=np.sum(bothinslens);basis['ALLdel_len_tot']=np.sum(bothdelslens)
    return basis


#%%
count_fixed_mutations_s1(r'C:\Users\rmoge\GitHub\MuriSam\MA_testing\3B_1_clean.VCF')
#%%
wholedict=dict()
for n in range(1,97):
    stringname="s1_" + str(n)
    try:
        thisfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_VCF\s1_' + str(n) + r'_clean.VCF'
        wholedict[stringname]=count_fixed_mutations_s1(thisfile)
    except FileNotFoundError:
        continue
MA_s1_df=pd.DataFrame.from_dict(wholedict,orient="index")
MA_s1_df.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_VCF\MA_s1_summary.csv') 

#%%
#Now, we need to count the number of nonsyn and syn SNMs in all of the lines.

def calc_dNdS_fromtuple(intup):
    N=intup[0];S=intup[1];DN=intup[2];DS=intup[3]
    dN = -3/4*math.log(1-(4/3*DN/N))
    dS = -3/4*math.log(1-(4/3*DS/S))
    #print("dN = " + str(dN)); print("dS = " + str(dS))
    #print("dN/dS = " + str(dN/dS))
    return dN/dS

def SNMcounter_s1(infile=""):
    fileN=0;fileS=1;unnamed=dict()#adding a pseudocount of S = 1 because some of the lines have S = 0, and I cannot divide by 0
    with open(infile, 'r') as ingd:
        for index,line in enumerate(csv.reader(ingd,delimiter='\t')):
 #           if index < 45 and index > 1:
            if line[0]!="SNP":
                continue
            elif line[4]=="917603" or line[4]=="31655" or line[4]=="120719" or line[4]=="359605" or line[4]=="549055" or line[4]=="118607" or line[4]=="538105" or line[4]=="855298":
                continue
            else:
                if line[6][0:2]=="aa" and line[13]=="frequency=1":
                    #print(str(line[6][-2:] + str(line[8][-2:])))
                    if line[6][-2:] == line[8][-2:]:
                        fileS+=1
                    else:
                        fileN+=1
                else:
                    continue
    #print(fileN);print(fileS)
    unnamed["dN"]=fileN;unnamed["dS"]=fileS
    
    return unnamed

'''
              if m.POS >= 76900 and m.POS <= 76909:#there is an oligo-T deletion in the MA_3B ancestor, but breseq doesn't always score it exactly the same, and therefore it wasn't subtracted out of the lines---hence I am subtracting it here
            continue
        elif m.POS == 20938:'''
def SNMcounter_3b(infile=""):
    fileN=0;fileS=1;unnamed=dict()#adding a pseudocount of S = 1 because some of the lines have S = 0, and I cannot divide by 0
    with open(infile, 'r') as ingd:
        for index,line in enumerate(csv.reader(ingd,delimiter='\t')):
 #           if index < 45 and index > 1:
            if line[0]!="SNP":
                continue
            elif (int(line[4]) >= 76900 and int(line[4]) <= 76909) or line[4]=="20938" or line[4]=="484801" or line[4]=="484793" or line[4]=="484795" or line[4]=="516354" or line[4]=="484791" or line[4]=="516356" or line[4]=="516362" or line[4]=="420163" or line[4]=="20023" or line[4]=="14906" or line[4]=="385785" or line[4]=="76905" or line[4]=="145853" or line[4]=="437910" or line[4]=="437909" or line[4]=="449149":
                continue
            else:
                #print(line[4])
                if line[6][0:2]=="aa" and line[13]=="frequency=1":
                    #print(str(line[6][-2:] + str(line[8][-2:])))
                    if line[6][-2:] == line[8][-2:]:
                        fileS+=1
                    else:
                        fileN+=1
                else:
                    continue
    #print(fileN);print(fileS)
    unnamed["dN"]=fileN;unnamed["dS"]=fileS
    
    return unnamed

def intergenic_or_coding_s1(infile=""):
    fileI=0;fileC=0;nameless=dict()
    with open(infile, 'r') as ingd:
        for index,line in enumerate(csv.reader(ingd,delimiter='\t')):
            if line[0]!="SNP":
                continue
            elif line[4]=="917603" or line[4]=="31655" or line[4]=="120719" or line[4]=="359605" or line[4]=="549055" or line[4]=="118607" or line[4]=="538105" or line[4]=="855298":
                continue
            else:
                if line[8][:11]=="aa_ref_seq=" and line[13]=="frequency=1":
                    #print('correct length')
                    fileC+=1
                elif line[8][:24]=="gene_position=intergenic" and line[6]=="frequency=1":
                    fileI+=1
                elif line[8][:23]=="gene_position=noncoding" and line[6]=="frequency=1":
                    fileI+=1
                else:
                    continue
          
    nameless['in_CDS']=fileC;nameless['not_in_CDS']=fileI
    return nameless
            
def intergenic_or_coding_3b(infile=""):
    fileI=0;fileC=0;nameless=dict()
    with open(infile, 'r') as ingd:
        for index,line in enumerate(csv.reader(ingd,delimiter='\t')):
            if line[0]!="SNP":
                continue
            elif (int(line[4]) >= 76900 and int(line[4]) <= 76909) or line[4]=="20938" or line[4]=="484801" or line[4]=="484793" or line[4]=="484795" or line[4]=="516354" or line[4]=="484791" or line[4]=="516356" or line[4]=="516362" or line[4]=="420163" or line[4]=="20023" or line[4]=="14906" or line[4]=="385785" or line[4]=="76905" or line[4]=="145853" or line[4]=="437910" or line[4]=="437909" or line[4]=="449149":
                continue
            else:
                if line[8][:11]=="aa_ref_seq=" and line[13]=="frequency=1":
                    #print('correct length')
                    fileC+=1
                elif line[8][:24]=="gene_position=intergenic" and line[6]=="frequency=1":
                    fileI+=1
                elif line[8][:23]=="gene_position=noncoding" and line[6]=="frequency=1":
                    fileI+=1
                else:
                    continue
          
    nameless['in_CDS']=fileC;nameless['not_in_CDS']=fileI
    return nameless            
    #are you gonna combine them all, and then do a big 2-proportion chisq test with continuity correction? Or test each line individually in the same way
 
#testfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_annotatedGD\s1_1_sub_annotated.gd'
#babycakes=intergenic_or_coding(testfile)
    

testfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_annotatedGD\s1_1_sub_annotated.gd'
babycakes=SNMcounter_3b(testfile)

#%% 
    
SNMdicts1=dict()
for n in range(1,97):
    stringname="s1_" + str(n)
    try:
        thisfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_annotatedGD\s1_' + str(n) + r'_sub_annotated.gd'
        SNMdicts1[stringname]=SNMcounter_s1(thisfile)
        inputter=tuple([2079386,544861,SNMdicts1[stringname]["dN"],SNMdicts1[stringname]["dS"]])
        #print(inputter)
        #print(SNMdicts1[stringname]["N"])
        SNMdicts1[stringname]["dN_dS"]=calc_dNdS_fromtuple(inputter)
        
        tempd=intergenic_or_coding_s1(thisfile)
        SNMdicts1[stringname]['in_CDS']=tempd['in_CDS'];SNMdicts1[stringname]['not_in_CDS']=tempd['not_in_CDS']
        
    except FileNotFoundError:
        continue
MA_s1_dnds=pd.DataFrame.from_dict(SNMdicts1,orient="index")
MA_s1_dnds.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_annotatedGD\MA_s1_dnds.csv')   
    
SNMdict3b=dict()
for n in range(1,97):
    stringname="3B_" + str(n)
    try:
        thisfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_3B_annotatedGD\3B_' + str(n) + r'_sub_annotated.gd'
        SNMdict3b[stringname]=SNMcounter_3b(thisfile)
        inputter2=tuple([1081335,287360,SNMdict3b[stringname]["dN"],SNMdict3b[stringname]["dS"]])
        SNMdict3b[stringname]["dN_dS"]=calc_dNdS_fromtuple(inputter2)
        
        tempd=intergenic_or_coding_3b(thisfile)
        SNMdict3b[stringname]['in_CDS']=tempd['in_CDS'];SNMdict3b[stringname]['not_in_CDS']=tempd['not_in_CDS']
        
    except FileNotFoundError:
        continue
MA_s1_dnds=pd.DataFrame.from_dict(SNMdict3b,orient="index")
MA_s1_dnds.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_3B_annotatedGD\MA_3B_dnds.csv')      
#%%

    
    
    
    
    
    
    
    
    
    
from __future__ import division
import os, math, pickle
import mult_tools as mt
import numpy as np
from scipy.stats import ttest_ind
import pandas as pd

import matplotlib.pyplot as plt
from decimal import Decimal
from matplotlib.ticker import FormatStrFormatter

mydir = os.path.expanduser("~/Desktop/OUTPUTS")
output_to_keep = ['INS', 'DEL', 'SNP', 'SUB']
strains = ['wildtype', 'minimal']


# two arguments: minimal or wildtype

def get_sites_to_remove(strain):
    if strain == 'minimal':
        acnestor_path = mydir + '/syn3B_minimal/mmW_3B.ancestor/output.gd'
    elif strain == 'wildtype':
        acnestor_path = mydir + '/syn1.0_wildtype/mm8_syn1.0.ancestor/output.gd'
    sites_to_remove = []
    for i, line in enumerate(open(acnestor_path, 'r')):
        line_split = line.strip().split('\t')
        if line_split[0] in output_to_keep:
            sites_to_remove.append( line_split[3] + '_' + str(line_split[4]))
    return sites_to_remove



def get_multiplicity(nmin = 2, FDR = 0.05):
    p_star_dict = {}
    G_score_list = []
    for strain in strains:

        sites_to_remove = get_sites_to_remove(strain)
        gene_count_dict = {}
        if strain == 'minimal':
            dirs = ['syn3B_minimal/mm13', 'syn3B_minimal/mm11', 'syn3B_minimal/mm10', 'syn3B_minimal/mm9']
            ref_path = mydir + '/syn3B_minimal/reference/Synthetic.bacterium_JCVI-Syn3A.gb'
        elif strain == 'wildtype':
            dirs = ['syn1.0_wildtype/mm6', 'syn1.0_wildtype/mm4', 'syn1.0_wildtype/mm3', 'syn1.0_wildtype/mm1']
            ref_path = mydir + '/syn1.0_wildtype/reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb'
        effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            muts = 0
            for i, line in enumerate(open(mydir+'/'+dir+'/annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                if line_split[0] == 'SNP':
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous':
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict:
                                    gene_count_dict[locus_tag_j] = 0
                                gene_count_dict[locus_tag_j] += 1
                        else:
                            if locus_tag not in gene_count_dict:
                                gene_count_dict[locus_tag] = 0
                            gene_count_dict[locus_tag] += 1

                    else:
                        continue
                else:
                    if len([s for s in line_split if 'gene_position=coding' in s]) >= 1:
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict:
                                    gene_count_dict[locus_tag_j] = 0
                                gene_count_dict[locus_tag_j] += 1

                        else:
                            if locus_tag not in gene_count_dict:
                                gene_count_dict[locus_tag] = 0
                            gene_count_dict[locus_tag] += 1

        # get multiplicity scores
        gene_parallelism_statistics = {}
        for gene_i, length_i in effective_gene_lengths.items():
            gene_parallelism_statistics[gene_i] = {}
            gene_parallelism_statistics[gene_i]['length'] = length_i
            gene_parallelism_statistics[gene_i]['observed'] = 0
            gene_parallelism_statistics[gene_i]['multiplicity'] = 0

        # save number of mutations for multiplicity
        for locus_tag_i, n_muts_i in gene_count_dict.items():
            gene_parallelism_statistics[locus_tag_i]['observed'] = n_muts_i

        L_mean = np.mean(list(effective_gene_lengths.values()))
        L_tot = sum(list(effective_gene_lengths.values()))
        n_tot = sum(gene_count_dict.values())
        # don't include taxa with less than 20 mutations
        print("N_total = " + str(n_tot))
        # go back over and calculate multiplicity
        for locus_tag_i in gene_parallelism_statistics.keys():
            # double check the measurements from this
            gene_parallelism_statistics[locus_tag_i]['multiplicity'] = gene_parallelism_statistics[locus_tag_i]['observed'] *1.0/ effective_gene_lengths[locus_tag_i] * L_mean
            gene_parallelism_statistics[locus_tag_i]['expected'] = n_tot*gene_parallelism_statistics[locus_tag_i]['length']/L_tot

        pooled_multiplicities = np.array([gene_parallelism_statistics[gene_name]['multiplicity'] for gene_name in gene_parallelism_statistics.keys() if gene_parallelism_statistics[gene_name]['multiplicity'] >=1])
        pooled_multiplicities.sort()

        pooled_tupe_multiplicities = np.array([(gene_parallelism_statistics[gene_name]['multiplicity'], gene_parallelism_statistics[gene_name]['observed']) for gene_name in gene_parallelism_statistics.keys() if gene_parallelism_statistics[gene_name]['multiplicity'] >=1])
        pooled_tupe_multiplicities = sorted(pooled_tupe_multiplicities, key=lambda x: x[0])
        pooled_tupe_multiplicities_x = [i[0] for i in pooled_tupe_multiplicities]
        pooled_tupe_multiplicities_y = [i[1] for i in pooled_tupe_multiplicities]
        pooled_tupe_multiplicities_y = [sum(pooled_tupe_multiplicities_y[i:]) / sum(pooled_tupe_multiplicities_y) for i in range(len(pooled_tupe_multiplicities_y))]

        null_multiplicity_survival = mt.NullGeneMultiplicitySurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )
        null_multiplicity_survival_copy = null_multiplicity_survival(pooled_multiplicities)
        null_multiplicity_survival_copy = [sum(null_multiplicity_survival_copy[i:]) / sum(null_multiplicity_survival_copy) for i in range(len(null_multiplicity_survival_copy)) ]
        #threshold_idx = numpy.nonzero((null_multiplicity_survival(observed_ms)*1.0/observed_multiplicity_survival)<FDR)[0][0]
        mult_survival_dict = {'Mult': pooled_multiplicities, 'Obs_fract': pooled_tupe_multiplicities_y, 'Null_fract': null_multiplicity_survival_copy}
        mult_survival_df = pd.DataFrame(mult_survival_dict)
        mult_survival_df_out = mydir + '/mult_survival_curves_' + strain + '.txt'
        mult_survival_df.to_csv(mult_survival_df_out, sep = '\t', index = True)


        # get likelihood score and null test
        observed_G, pvalue = mt.calculate_total_parallelism(gene_parallelism_statistics)
        G_score_list.append((strain, observed_G, pvalue))
        print(strain, observed_G, pvalue)

        # Give each gene a p-value, get distribution
        gene_logpvalues = mt.calculate_parallelism_logpvalues(gene_parallelism_statistics)
        pooled_pvalues = []
        for gene_name in gene_logpvalues.keys():
            if (gene_parallelism_statistics[gene_name]['observed']>= nmin) and (float(gene_logpvalues[gene_name]) >= 0):
                pooled_pvalues.append( gene_logpvalues[gene_name] )

        pooled_pvalues = np.array(pooled_pvalues)
        pooled_pvalues.sort()
        if len(pooled_pvalues) == 0:
            continue


        null_pvalue_survival = mt.NullGeneLogpSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics, nmin=nmin)
        observed_ps, observed_pvalue_survival = mt.calculate_unnormalized_survival_from_vector(pooled_pvalues, min_x=-4)
        # Pvalue version
        # remove negative minus log p values.
        neg_p_idx = np.where(observed_ps>=0)
        observed_ps_copy = observed_ps[neg_p_idx]
        observed_pvalue_survival_copy = observed_pvalue_survival[neg_p_idx]
        pvalue_pass_threshold = np.nonzero(null_pvalue_survival(observed_ps_copy)*1.0/observed_pvalue_survival_copy<FDR)[0]
        if len(pvalue_pass_threshold) == 0:
            continue
        threshold_idx = pvalue_pass_threshold[0]
        pstar = observed_ps_copy[threshold_idx] # lowest value where this is true
        num_significant = observed_pvalue_survival[threshold_idx]
        # make it log base 10
        logpvalues_dict = {'P_value': observed_ps/math.log(10), 'Obs_num': observed_pvalue_survival, 'Null_num': null_pvalue_survival(observed_ps)}
        logpvalues_df = pd.DataFrame(logpvalues_dict)
        logpvalues_df_out = mydir + '/logpvalues_' + strain + '.txt'
        logpvalues_df.to_csv(logpvalues_df_out, sep = '\t', index = True)

        p_star_dict[strain] = (num_significant, pstar/math.log(10))


        output_mult_gene_filename = mydir + '/mult_genes_sig_' + strain + '.txt'
        output_mult_gene = open(output_mult_gene_filename,"w")
        output_mult_gene.write(",".join(["Gene", "Length", "Observed", "Expected", "Multiplicity", "-log10(P)"]))
        for gene_name in sorted(gene_parallelism_statistics, key=lambda x: gene_parallelism_statistics.get(x)['observed'],reverse=True):
            if gene_logpvalues[gene_name] >= pstar and gene_parallelism_statistics[gene_name]['observed']>=nmin:
                output_mult_gene.write("\n")
                # log base 10 transform the p-values here as well
                output_mult_gene.write("%s, %0.1f, %d, %0.2f, %0.2f, %g" % (gene_name, gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], abs(gene_logpvalues[gene_name])/math.log(10) ))
        output_mult_gene.close()


    total_parallelism_path = mydir + '/total_parallelism.txt'
    total_parallelism = open(total_parallelism_path,"w")
    total_parallelism.write("\t".join(["Strain", "G_score", "p_value"]))
    for i in range(len(G_score_list)):
        taxon_i = G_score_list[i][0]
        G_score_i = G_score_list[i][1]
        p_value_i = G_score_list[i][2]
        total_parallelism.write("\n")
        total_parallelism.write("\t".join([taxon_i, str(G_score_i), str(p_value_i)]))

    total_parallelism.close()
    with open(mydir + '/p_star.txt', 'wb') as file:
        file.write(pickle.dumps(p_star_dict)) # use `pickle.loads` to do the reverse



def plot_multiplicity_survival():
    df_par = pd.read_csv(mydir + '/total_parallelism.txt', sep = '\t' )
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(strains)):
        strain = strains[i]
        df_path = mydir + '/mult_survival_curves_' + strain + '.txt'
        df = pd.read_csv(df_path, sep = '\t', index_col=0)
        new_x = [1.0] + df.Mult.tolist() + [df.Mult.tolist()[-1]]
        new_obs_y =[1.0] + df.Obs_fract.tolist() + [ 0.0001]
        new_null_y = [1.0] + df.Null_fract.tolist() + [ 0.0001]

        ax = fig.add_subplot(2, 1, i+1)
        ax.plot(new_x, new_obs_y, '-', c='royalblue', lw=4, alpha = 0.8, zorder=1)
        ax.plot(new_x, new_null_y, '-', c='dimgrey', lw=4, alpha = 0.8, zorder=0)
        ax.set_xlim([0.9, 16])

        taxon_par = df_par.loc[df_par['Strain'] == strain]

        ax.annotate(r'$\Delta \ell= $'+ str(round(float(taxon_par.G_score), 3)), (0.6 *max(new_x), 0.9), fontsize=8)
        if np.log10(float(taxon_par.p_value)) < -3:
            ax.annotate(r'$\mathrm{p} = $'+ str('%.2E' % Decimal(float(taxon_par.p_value))), (0.6 *max(new_x), 0.75), fontsize=8)
        else:
            ax.annotate(r'$\mathrm{p} = $'+ str(round(float(taxon_par.p_value),3)), (0.6 *max(new_x), 0.75), fontsize=8)

        ax.title.set_text(strain)
        ax.title.set_fontsize(8.5)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    fig.text(0.5, 0.02, 'Gene multiplicity, ' + '$m$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Fraction mutations ' + '$\geq m$', va='center', rotation='vertical', fontsize=16)

    fig_name = mydir + '/mult_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def plot_logpvalue_survival():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    pstar_dict = pickle.load(open(mydir + '/p_star.txt', 'rb'))
    for i in range(0, len(strains)):
        strain = strains[i]
        pstar_i = pstar_dict[strain][1]
        num_significant_i = pstar_dict[strain][0] -1
        df = pd.read_csv(mydir + '/logpvalues_' + strain + '.txt', sep = '\t', index_col=0)
        new_x = df.P_value.tolist()
        new_obs_y = df.Obs_num.tolist()
        new_null_y = df.Null_num.tolist()

        ax = fig.add_subplot(2, 1, i+1)

        ax.plot(new_x, new_null_y, '-', c='dimgrey', lw=4, alpha = 0.8, zorder=0)
        ax.plot(new_x, new_obs_y, '-', c='royalblue', lw=4, alpha = 0.8, zorder=1)
        if pstar_i <0:
            y_range = [f[1] for f in list(zip(new_x, new_obs_y)) if f[0] > 0]
            ax.plot([1, 1],[5e-02,max(y_range)],'k-',linewidth=0.5, zorder=2)
            ax.plot([-3,1],[max(y_range), max(y_range)],'k-',linewidth=0.5, zorder=3)
            ax.plot([1], [max(y_range)], c='r', marker='o', zorder=4)
        else:
            ax.plot([pstar_i, pstar_i],[5e-02,num_significant_i],'k-',linewidth=0.5, zorder=2)
            ax.plot([-3,pstar_i],[num_significant_i, num_significant_i],'k-',linewidth=0.5, zorder=3)
            ax.plot([pstar_i], [num_significant_i], c='r', marker='o', zorder=4)

        ax.set_xlim([0.25, 13])

        ax.title.set_text(strain)

        ax.title.set_fontsize(8.5)

    fig.text(0.5, 0.02, '$-\mathrm{log}_{10}P$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Number of genes', va='center', rotation='vertical', fontsize=16)

    fig_name = mydir + '/logpvalue_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




# equal rate of evolution
#print(ttest_ind([14,23,15,20], [19,14,13,9],equal_var=False))


#get_multiplicity('wildtype')
#plot_multiplicity_survival()
plot_logpvalue_survival()

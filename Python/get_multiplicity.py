from __future__ import division
import os, math, pickle
import mult_tools as mt
import numpy as np
import pandas as pd

from decimal import Decimal

output_to_keep = ['INS', 'DEL', 'SNP', 'SUB']
strains = ['wildtype', 'minimal']


# two arguments: minimal or wildtype

def get_sites_to_remove(strain):
    if strain == 'minimal':
        acnestor_path = mt.get_path() + '/data/syn3B_minimal/mmW_3B.ancestor/output.gd'
    elif strain == 'wildtype':
        acnestor_path = mt.get_path() + '/data/syn1.0_wildtype/mm8_syn1.0.ancestor/output.gd'
    sites_to_remove = []
    for i, line in enumerate(open(acnestor_path, 'r')):
        line_split = line.strip().split('\t')
        if line_split[0] in output_to_keep:
            sites_to_remove.append( line_split[3] + '_' + str(line_split[4]))
    return sites_to_remove



def get_multiplicity(nmin = 2, FDR = 0.05):
    p_star_dict = {}
    G_score_list = []

    gene_by_pop_dict = {}
    for strain in strains:

        sites_to_remove = get_sites_to_remove(strain)
        gene_count_dict = {}
        if strain == 'minimal':
            dirs = ['syn3B_minimal/mm13', 'syn3B_minimal/mm11', 'syn3B_minimal/mm10', 'syn3B_minimal/mm9']
            ref_path = mt.get_path() + '/data/syn3B_minimal/reference/Synthetic.bacterium_JCVI-Syn3A.gb'
        elif strain == 'wildtype':
            dirs = ['syn1.0_wildtype/mm6', 'syn1.0_wildtype/mm4', 'syn1.0_wildtype/mm3', 'syn1.0_wildtype/mm1']
            ref_path = mt.get_path() + '/data/syn1.0_wildtype/reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb'
        effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            for i, line in enumerate(open(mt.get_path()+'/data/'+dir+'/annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                if frequency != 1:
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
        mult_survival_df_out = mt.get_path() + '/data/mult_survival_curves_' + strain + '.txt'
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
        logpvalues_df_out = mt.get_path() + '/data/logpvalues_' + strain + '.txt'
        logpvalues_df.to_csv(logpvalues_df_out, sep = '\t', index = True)

        p_star_dict[strain] = (num_significant, pstar/math.log(10))


        output_mult_gene_filename = mt.get_path() + '/data/mult_genes_sig_' + strain + '.txt'
        output_mult_gene = open(output_mult_gene_filename,"w")
        output_mult_gene.write(",".join(["Gene", "Length", "Observed", "Expected", "Multiplicity", "-log10(P)"]))
        for gene_name in sorted(gene_parallelism_statistics, key=lambda x: gene_parallelism_statistics.get(x)['observed'],reverse=True):
            if gene_logpvalues[gene_name] >= pstar and gene_parallelism_statistics[gene_name]['observed']>=nmin:
                output_mult_gene.write("\n")
                # log base 10 transform the p-values here as well
                output_mult_gene.write("%s, %0.1f, %d, %0.2f, %0.2f, %g" % (gene_name, gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], abs(gene_logpvalues[gene_name])/math.log(10) ))
        output_mult_gene.close()


    total_parallelism_path = mt.get_path() + '/data/total_parallelism.txt'
    total_parallelism = open(total_parallelism_path,"w")
    total_parallelism.write("\t".join(["Strain", "G_score", "p_value"]))
    for i in range(len(G_score_list)):
        taxon_i = G_score_list[i][0]
        G_score_i = G_score_list[i][1]
        p_value_i = G_score_list[i][2]
        total_parallelism.write("\n")
        total_parallelism.write("\t".join([taxon_i, str(G_score_i), str(p_value_i)]))

    total_parallelism.close()
    with open(mt.get_path() + '/data/p_star.txt', 'wb') as file:
        file.write(pickle.dumps(p_star_dict)) # use `pickle.loads` to do the reverse


def get_multiplicity_matrix():

    gene_by_pop_dict = {}#initialize GxP matrix
    for strain in strains:
        sites_to_remove = get_sites_to_remove(strain)
        if strain == 'minimal':
            dirs = ['syn3B_minimal/mm13', 'syn3B_minimal/mm11', 'syn3B_minimal/mm10', 'syn3B_minimal/mm9']
            ref_path = mt.get_path() + '/data/syn3B_minimal/reference/Synthetic.bacterium_JCVI-Syn3A.gb'
        elif strain == 'wildtype':
            dirs = ['syn1.0_wildtype/mm6', 'syn1.0_wildtype/mm4', 'syn1.0_wildtype/mm3', 'syn1.0_wildtype/mm1']
            ref_path = mt.get_path() + '/data/syn1.0_wildtype/reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb'

        effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            pop = dir.split('/')[1]
            gene_count_dict_pop = {}
            gene_by_pop_dict[pop] = {}#GxP matrix add "rows", i.e., sub-dicts
            for i, line in enumerate(open(mt.get_path()+'/data/'+dir+'/annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                if frequency != 1:
                    continue
                if line_split[0] == 'SNP':
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous':
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += 1
                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += 1

                    else:
                        continue
                else:
                    if len([s for s in line_split if 'gene_position=coding' in s]) >= 1:
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += 1

                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += 1

            gene_parallelism_statistics = {}
            for gene_i, length_i in effective_gene_lengths.items():
                gene_parallelism_statistics[gene_i] = {}
                gene_parallelism_statistics[gene_i]['length'] = length_i
                gene_parallelism_statistics[gene_i]['observed'] = 0
                gene_parallelism_statistics[gene_i]['multiplicity'] = 0

            # save number of mutations for multiplicity
            for locus_tag_i, n_muts_i in gene_count_dict_pop.items():
                gene_parallelism_statistics[locus_tag_i]['observed'] = n_muts_i

            # save number of mutations for multiplicity
            L_mean = np.mean(list(effective_gene_lengths.values()))
            L_tot = sum(list(effective_gene_lengths.values()))
            n_tot = sum(gene_count_dict_pop.values())
            # go back over and calculate multiplicity
            for locus_tag_i in gene_parallelism_statistics.keys():
                # double check the measurements from this
                gene_parallelism_statistics[locus_tag_i]['multiplicity'] = gene_parallelism_statistics[locus_tag_i]['observed'] *1.0/ effective_gene_lengths[locus_tag_i] * L_mean
                gene_parallelism_statistics[locus_tag_i]['expected'] = n_tot*gene_parallelism_statistics[locus_tag_i]['length']/L_tot

            # split locus tags
            for locus_tag_i in gene_parallelism_statistics.keys():
                mult_i = gene_parallelism_statistics[locus_tag_i]['multiplicity']
                if mult_i > 0:
                    locus_tag_i_num = locus_tag_i.split('_')[1]
                    gene_by_pop_dict[pop][locus_tag_i_num] = mult_i

    gene_by_pop_df = pd.DataFrame(gene_by_pop_dict)#outputs the GxP matrix
    gene_by_pop_df = gene_by_pop_df.T
    gene_by_pop_df.fillna(0, inplace=True)

    gene_by_pop_df_out = mt.get_path() + '/data/mult_by_pop.txt'
    gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = '\t', index = True)




#get_multiplicity_matrix()
# equal rate of evolution

#get_multiplicity()

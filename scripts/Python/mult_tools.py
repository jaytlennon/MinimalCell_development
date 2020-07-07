from __future__ import division
import os
from Bio import SeqIO
import numpy as np
from scipy.stats import poisson
from scipy.special import gammaln
from scipy.linalg import block_diag

#from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances

def get_path():
    return os.path.expanduser("~/GitHub/MinimalCell")

def get_pop_dict():
    # min mininmal wt wildtype
    pop_dict = {'mm13':'min', 'mm11':'min', 'mm10':'min', 'mm9':'min',
                'mm6':'wt', 'mm4':'wt', 'mm3':'wt','mm1':'wt'}
    return pop_dict




def get_F_2(X, N1, N2):
    '''
    Modified F-statistic from Anderson et al., 2017 doi: 10.1111/anzs.12176
    Function assumes that the rows of the count matrix are sorted by group
    i.e., group one is first N1 rows, rest of the N2 rows are group two
    '''
    N = N1+N2
    dist_matrix = euclidean_distances(X, X)
    A = -(1/2) * (dist_matrix ** 2)
    I = np.identity(N)
    J_N = np.full((N, N), 1)
    G = (I - ((1/N) * J_N )) @ A @ (I - ((1/N) * J_N ))
    n1 = (1/N1) * np.full((N1, N1), 1)
    n2 = (1/N2) * np.full((N2, N2), 1)
    H = block_diag(n1, n2) - ((1/N) * J_N )
    # indicator matrices
    U_1 = np.diag( (N1*[1]) + (N2*[0]))
    U_2 = np.diag( (N1*[0]) + (N2*[1]))

    V_1 = np.trace(((I - H) @ U_1 @ (I - H)) @ G ) / (N1-1)
    V_2 = np.trace(((I - H) @ U_2 @ (I - H)) @ G ) / (N2-1)

    F_2 = np.trace(H @ G) / (((1- (N1/N)) * V_1) + ((1- (N2/N)) * V_2))

    return F_2, V_1, V_2


# NullMultiplicitySurvivalFunction class is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
class NullGeneMultiplicitySurvivalFunction(object):
    # Null multiplicity distribution for genes

    def __init__(self, Ls, ntot):
        self.ntot = ntot
        self.Ls = np.array(Ls)
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls*1.0/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps

    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics):
        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']

        return cls(Ls, ntot)

    def __call__(self, m):
        #lower_limits = np.ceil(m[:,None]*self.Ls[None,:]/self.Lavg)-1+0.1
        #return (poisson.sf(lower_limits, self.expected_ns[None,:])).sum(axis=1)
        lower_limits = np.ceil(m[:,None]*self.Ls[None,:]/self.Lavg)-2+0.1
        return (poisson.sf(lower_limits, self.expected_ns[None,:])*self.ps[None,:]).sum(axis=1)





# calculate_unnormalized_survival_from_vector function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_unnormalized_survival_from_vector(xs, min_x=None, max_x=None, min_p=1e-10):
    if min_x==None:
        min_x = xs.min()-1

    if max_x==None:
        max_x = xs.max()+1

    unique_xs = set(xs)
    unique_xs.add(min_x)
    unique_xs.add(max_x)

    xvalues = []
    num_observations = []

    for x in sorted(unique_xs):
        xvalues.append(x)
        num_observations.append( (xs>=x).sum() )

    # So that we can plot CDF, SF on log scale
    num_observations[0] -= min_p
    num_observations[1] -= min_p
    num_observations[-1] += min_p

    return np.array(xvalues), np.array(num_observations)

# calculate_G_scores function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_G_scores(gene_statistics, allowed_genes=None):
    # Calculates the G score for the whole gene, i.e.
    # n*g

    gene_g_scores = calculate_g_scores(gene_statistics,allowed_genes)

    gene_G_scores = {gene_name: gene_statistics[gene_name]['observed']*gene_g_scores[gene_name] for gene_name in gene_g_scores.keys()}

    return gene_G_scores


# calculate_total_parallelism function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_total_parallelism(gene_statistics, allowed_genes=None, num_bootstraps=10000):

    if allowed_genes==None:
        allowed_genes = gene_statistics.keys()

    Ls = []
    ns = []

    for gene_name in allowed_genes:

        Ls.append( gene_statistics[gene_name]['length'] )
        ns.append( gene_statistics[gene_name]['observed'] )


    Ls = np.array(Ls)
    ns = np.array(ns)

    Ltot = Ls.sum()
    ntot = ns.sum()
    ps = Ls*1.0/Ltot

    gs = ns*np.log(ns/(ntot*ps)+(ns==0))

    observed_G = gs.sum()/ns.sum()
    bootstrapped_Gs = []
    for bootstrap_idx in range(0,num_bootstraps):
        bootstrapped_ns = np.random.multinomial(ntot,ps)
        bootstrapped_gs = bootstrapped_ns*np.log(bootstrapped_ns/(ntot*ps)+(bootstrapped_ns==0))
        bootstrapped_G = bootstrapped_gs.sum()/bootstrapped_ns.sum()

        bootstrapped_Gs.append(bootstrapped_G)

    bootstrapped_Gs = np.array(bootstrapped_Gs)

    pvalue = ((bootstrapped_Gs>=observed_G).sum()+1.0)/(len(bootstrapped_Gs)+1.0)
    return observed_G, pvalue


# calculate_poisson_log_survival function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_poisson_log_survival(ns, expected_ns):

    survivals = poisson.sf(ns-0.1, expected_ns)

    logsurvivals = np.zeros_like(survivals)
    logsurvivals[survivals>1e-20] = -np.log(survivals[survivals>1e-20])
    logsurvivals[survivals<=1e-20] = (-ns*np.log(ns/expected_ns+(ns==0))+ns-expected_ns)[survivals<=1e-20]

    return logsurvivals


# calculate_parallelism_logpvalues function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_parallelism_logpvalues(gene_statistics):

    gene_names = []
    Ls = []
    ns = []
    expected_ns = []

    for gene_name in gene_statistics.keys():
        gene_names.append(gene_name)
        ns.append(gene_statistics[gene_name]['observed'])
        expected_ns.append(gene_statistics[gene_name]['expected'])

    ns = np.array(ns)
    expected_ns = np.array(expected_ns)

    logpvalues = calculate_poisson_log_survival(ns, expected_ns)

    return {gene_name: logp for gene_name, logp in zip(gene_names, logpvalues)}


# NullGeneLogpSurvivalFunction class is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
class NullGeneLogpSurvivalFunction(object):
    # Null distribution of -log p for each gene

    def __init__(self, Ls, ntot,nmin=0):
        self.ntot = ntot
        self.Ls = np.array(Ls)*1.0
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps
        self.nmin = nmin

    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics,nmin=0):

        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']

        return cls(Ls, ntot, nmin)

    def __call__(self, mlogps):

        # Do sum by hand
        ns = np.arange(0,400)*1.0

        logpvalues = calculate_poisson_log_survival(ns[None,:], self.expected_ns[:,None])

        logprobabilities = ns[None,:]*np.log(self.expected_ns)[:,None]-gammaln(ns+1)[None,:]-self.expected_ns[:,None]
        probabilities = np.exp(logprobabilities)
        survivals = np.array([ ((logpvalues>=mlogp)*(ns[None,:]>=self.nmin)*probabilities).sum() for mlogp in mlogps])
        return survivals





# calculate_synonymous_nonsynonymous_target_sizes function is modified from GitHub repo
# benjaminhgood/LTEE-metagenomic under GPL v2
def calculate_synonymous_nonsynonymous_target_sizes(gbf_path):
    codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R',
                    'CGC': 'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
                    'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C',
                    'TGC':'D', 'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E',
                    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H',
                    'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L',
                    'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                    'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F',
                    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'TCT':'S',
                    'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
                    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TGG':'W',
                    'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V',
                    'GTG':'V', 'TAA':'!', 'TGA':'!', 'TAG':'!' }

    # calculate number of synonymous opportunities for each codon
    codon_synonymous_opportunity_table = {}
    for codon in codon_table.keys():
        codon_synonymous_opportunity_table[codon] = {}
        for i in range(0,3):
            codon_synonymous_opportunity_table[codon][i] = -1 # since G->G is by definition synonymous, but we don't want to count it
            codon_list = list(codon)
            for base in ['A','C','T','G']:
                codon_list[i]=base
                new_codon = "".join(codon_list)
                if codon_table[codon]==codon_table[new_codon]:
                    # synonymous!
                    codon_synonymous_opportunity_table[codon][i]+=1

    codon_synonymous_substitution_table = {}
    codon_nonsynonymous_substitution_table = {}
    for codon in codon_table.keys():
        codon_synonymous_substitution_table[codon] = [[],[],[]]
        codon_nonsynonymous_substitution_table[codon] = [[],[],[]]

        for i in range(0,3):
            reference_base = codon[i]

            codon_list = list(codon)
            for derived_base in ['A','C','T','G']:
                if derived_base==reference_base:
                    continue
                substitution = '%s->%s' % (reference_base, derived_base)
                codon_list[i]=derived_base
                new_codon = "".join(codon_list)
                if codon_table[codon]==codon_table[new_codon]:
                    # synonymous!
                    codon_synonymous_substitution_table[codon][i].append(substitution)
                else:
                    codon_nonsynonymous_substitution_table[codon][i].append(substitution)
    bases = set(['A','C','T','G'])
    substitutions = []
    for b1 in bases:
        for b2 in bases:
            if b2==b1:
                continue
            substitutions.append( '%s->%s' % (b1,b2) )
    substitution_specific_synonymous_sites = {substitution: 0 for substitution in substitutions}
    substitution_specific_nonsynonymous_sites = {substitution: 0 for substitution in substitutions}
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}
    gene_length_map = {}

    if 'JCVI-Syn3A' in gbf_path:
        genome_size = 543379
    elif 'JCVI-syn1.0_CP002027' in gbf_path:
        genome_size = 1078809
    for record in SeqIO.parse(gbf_path, "genbank"):
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            if 'note' in feature.qualifiers:
                if 'incomplete' in feature.qualifiers['note'][0]:
                    continue
                if 'frameshifted' in feature.qualifiers['note'][0]:
                    continue
                if 'internal stop' in feature.qualifiers['note'][0]:
                    continue
            gene_name = feature.qualifiers['locus_tag'][0]

            if gene_name not in effective_gene_synonymous_sites:
                effective_gene_synonymous_sites[gene_name]=0
                effective_gene_nonsynonymous_sites[gene_name]=0
            # if there's no translation, ignore
            if 'translation' not in feature.qualifiers:
                continue
            aa_str = str(feature.qualifiers['translation'][0])
            nuc_str = str(feature.location.extract(record).seq[:-3])
            gene_length_map[gene_name] = len(nuc_str)
            for position in range(len(nuc_str)):
                codon_start = int(position/3)*3
                codon = nuc_str[codon_start:codon_start+3]
                if len(codon) <3:
                    continue
                position_in_codon = position%3

                effective_gene_synonymous_sites[gene_name] += codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                effective_gene_nonsynonymous_sites[gene_name] += 1-codon_synonymous_opportunity_table[codon][position_in_codon]/3.0

                for substitution in codon_synonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_synonymous_sites[substitution] += 1

                for substitution in codon_nonsynonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_nonsynonymous_sites[substitution] += 1

    substitution_specific_synonymous_fraction = {substitution: substitution_specific_synonymous_sites[substitution]*1.0/(substitution_specific_synonymous_sites[substitution]+substitution_specific_nonsynonymous_sites[substitution]) for substitution in substitution_specific_synonymous_sites.keys()}
    effective_gene_lengths = {gene_name: gene_length_map[gene_name]-effective_gene_synonymous_sites[gene_name] for gene_name in gene_length_map.keys()}
    effective_gene_lengths_synonymous = sum([effective_gene_synonymous_sites[gene_name] for gene_name in gene_length_map.keys()])
    effective_gene_lengths_nonsynonymous = sum([effective_gene_nonsynonymous_sites[gene_name] for gene_name in gene_length_map.keys()])
    effective_gene_lengths_noncoding = genome_size - effective_gene_lengths_synonymous-effective_gene_lengths_nonsynonymous
    # first two objects returned are dicts
    return effective_gene_lengths, effective_gene_synonymous_sites, effective_gene_lengths_synonymous, effective_gene_lengths_nonsynonymous, effective_gene_lengths_noncoding

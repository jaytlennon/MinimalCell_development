from __future__ import division
import pickle
import mult_tools as mt
import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

from sklearn.metrics import pairwise_distances
from skbio.stats.ordination import pcoa



def plot_multiplicity_survival():
    df_par = pd.read_csv(mt.get_path() + '/data/total_parallelism.txt', sep = '\t' )
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    for i in range(0, len(strains)):
        strain = strains[i]
        df_path = mt.get_path() + '/data/mult_survival_curves_' + strain + '.txt'
        df = pd.read_csv(df_path, sep = '\t', index_col=0)
        new_x = [1.0] + df.Mult.tolist() + [df.Mult.tolist()[-1]]
        new_obs_y =[1.0] + df.Obs_fract.tolist() + [ 0.0001]
        new_null_y = [1.0] + df.Null_fract.tolist() + [ 0.0001]

        ax = fig.add_subplot(2, 1, i+1)
        ax.plot(new_x, new_obs_y, '-', c='royalblue', lw=4, alpha = 0.8, zorder=1)
        ax.plot(new_x, new_null_y, '-', c='dimgrey', lw=4, alpha = 0.8, zorder=0)
        ax.set_xlim([0.9, 9])

        taxon_par = df_par.loc[df_par['Strain'] == strain]

        ax.annotate(r'$\Delta \ell= $'+ str(round(float(taxon_par.G_score), 3)), (0.6 *9, 0.9), fontsize=8)
        if np.log10(float(taxon_par.p_value)) < -3:
            ax.annotate(r'$\mathrm{p} = $'+ str('%.2E' % Decimal(float(taxon_par.p_value))), (0.6 *9, 0.75), fontsize=8)
        else:
            ax.annotate(r'$\mathrm{p} = $'+ str(round(float(taxon_par.p_value),3)), (0.6 *9, 0.75), fontsize=8)

        if strain == 'wildtype':
            ax.title.set_text('Wildtype')
        elif strain == 'minimal':
            ax.title.set_text('Minimal')
        ax.title.set_fontsize(12)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    fig.text(0.5, 0.02, 'Gene multiplicity, ' + '$m$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Fraction mutations ' + '$\geq m$', va='center', rotation='vertical', fontsize=16)

    fig_name = mt.get_path() + '/figures/mult_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def plot_logpvalue_survival():
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    pstar_dict = pickle.load(open(mt.get_path() + '/data/p_star.txt', 'rb'))
    for i in range(0, len(strains)):
        strain = strains[i]
        pstar_i = pstar_dict[strain][1]
        num_significant_i = pstar_dict[strain][0] -1
        df = pd.read_csv(mt.get_path() + '/data/logpvalues_' + strain + '.txt', sep = '\t', index_col=0)
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

        ax.set_xlim([0.25, 8])

        ax.title.set_text(strain)
        ax.title.set_fontsize(12)

    fig.text(0.5, 0.02, '$-\mathrm{log}_{10}P$', ha='center', fontsize=16)
    fig.text(0.02, 0.5, 'Number of genes', va='center', rotation='vertical', fontsize=16)

    fig_name = mt.get_path() + '/figures/logpvalue_survival.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)



def plot_pcoa(bs_iter = 10000):
    # get F stat and p value
    df = pd.read_csv(mt.get_path() + '/data/mult_by_pop.txt', sep = '\t', index_col=0)
    #mt.get_F_2(df,4,4)
    df = df/df.sum(axis=1)[:,None]
    df_bc = pairwise_distances(df, metric='braycurtis')

    df_pcoa = pcoa(df_bc , number_of_dimensions=3)
    ord_matrix = df_pcoa.samples

    F = mt.get_F_2(ord_matrix, 4,4)
    F_nulls = []
    for i in range(bs_iter):
        F_nulls.append(mt.get_F_2(ord_matrix.sample(frac=1), 4,4)[0])
    p_value = len([F_null for F_null in F_nulls if  F_null > F[0]]) / bs_iter
    print("F = " + str(round(F[0], 4)))
    print("p = " + str(round(p_value, 4)))

    #fig = plt.figure()
    fig, ax = plt.subplots(figsize=(6, 6))
    # Scatterplot on main ax
    ax.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
    ax.axvline(x=0, color='k', linestyle=':', alpha = 0.8, zorder=2)
    ax.scatter(0, 0, marker = "o", edgecolors='none', c = 'darkgray', s = 120, zorder=3)
    ax.scatter(ord_matrix.ix[0:4,0],ord_matrix.ix[0:4,1], marker = "o",
        edgecolors='#244162', c = 'blue', alpha = 0.8, s = 120, zorder=4, label='Wildtype')

    ax.scatter(ord_matrix.ix[4:,0],ord_matrix.ix[4:,1], marker = "o",
        edgecolors='#244162', c = 'r', alpha = 0.8, s = 120, zorder=4, label='Minimal cell')


    confidence_ellipse(ord_matrix.ix[0:4,0],ord_matrix.ix[0:4,1], ax,
        n_std=2, edgecolor='blue', linestyle='--', lw=3)
    confidence_ellipse(ord_matrix.ix[4:,0],ord_matrix.ix[4:,1], ax,
        n_std=2, edgecolor='red', linestyle='--', lw=3)
    #ax1.xlim([-0.7,0.7])
    #ax1.set_ylim([-0.7,0.7])

    ax.set_xlabel('PCo 1 (' + str(round(df_pcoa.proportion_explained[0],3)*100) + '%)' , fontsize = 14)
    ax.set_ylabel('PCo 2 (' + str(round(df_pcoa.proportion_explained[1],3)*100) + '%)' , fontsize = 14)



    plt.legend(loc="upper right")

    fig_name = mt.get_path() + '/figures/pcoa.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


plot_pcoa()

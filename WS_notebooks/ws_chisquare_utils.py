# ws_chisquare_utils.py
# WESmith 12/07/22
# tailored utility functions for chi-square experimentation

# Copyright (c) 2022 Warren E Smith  smiwarsky@gmail.com
# MIT License


from scipy.stats import chi2
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import pandas as pd


def plot_chi_square(n=6, xlim=10, ylim=1, fx=10, fy=6):
    fig, ax = plt.subplots(1, 1, figsize=(fx, fy))
    x = np.linspace(0, xlim, 1000)
    for df in range(1, n+1):
        ax.plot(x, chi2.pdf(x, df), label='df: {}'.format(df))
    ax.set_ylim(0.0, ylim)
    ax.grid()
    ax.legend()
    plt.show()


def plot_monte_carlo_chi_square(num=10000, bins=50, 
                                nr=5, nc=5, fx=15, fy=10):
    rows = list(range(nr))
    cols = list(range(nc))
    vals = [(r, c) for r in rows for c in cols]
    nn   = np.array(range(len(rows)*len(cols))) + 1
    fig, ax = plt.subplots(len(rows), len(cols), 
                           figsize=(fx, fy))
    for (r,c), n in zip(vals, nn):
        out = []
        for j in range(num):
            dd = np.random.normal(loc=0, scale=1.0, size=n)
            # form sum of squares of n normal RVs
            out.append(np.sum(dd*dd))
        out = np.array(out)
        ax[r, c].hist(out, bins=bins, density=True)
        x = np.linspace(out.min(), out.max(), 1000)
        ax[r, c].plot(x, chi2.pdf(x, n), 'r', 
                      label='df: {}'.format(n))
        if n == 1: ax[r, c].set_ylim(0.0, 1.0)
        ax[r, c].legend()
        ax[r, c].grid()
    
    txt = ('\nCHI-SQUARE DISTRIBUTION FOR {} '
           'SAMPLES OVER {} BINS\n'.format(num, bins))
    fig.suptitle(txt, fontsize='xx-large')
    plt.tight_layout()
    plt.show()


def create_chi_square_table(rows=None, cols=None, dataframe=False):
    '''
    Match chi-square table in Schaum's 'Statistics' p.345
    Verified with spot checks 12/7/22
    rows: list of desired degrees of freedom (default is Schaum's values)
    cols: list of desired percentile values (default is Schaum's values)
          eg, a percentile value of 0.95 returns the chi-sqare
           value that has 95% of the area to its LEFT; subsequently this
           same value represents 5% of the area to its RIGHT
    dataframe: BOOL if True, return a pandas dataframe, 
               otherwise (default) return
               a double dictionary indexed as follows:
               [degrees_of_freedom][percentile_value]: 
               eg, out[5][.95] would return 11.070497693516351
    '''
    default_rows = list(range(1, 101))
    default_cols = [.995, .99, .975, .95, .90, .75, 
                    .50, .25, .10, .05, .025, .01, .005]
    rows = rows if rows else default_rows
    cols = cols if cols else default_cols
    out = defaultdict()
    for df in rows:
        tmp = defaultdict()
        for percent in cols:
            tmp[percent] = chi2.ppf(percent, df)
        out[df] = tmp
    return pd.DataFrame(out).T if dataframe else out


def get_pvals(n, random=0):
    # create test p values for a multinomial distribution;
    # random = 0 returns a uniform distribution (ie, a fair die)
    p = np.random.normal(loc=0, scale=random, size=n) + 1.0
    if np.any(p <= 0):
        txt = ('getting negative or zero p values: '
               'make "random" smaller than {}'.format(random))
        raise ValueError(txt)
    return p / np.sum(p)


def chi_square_demo(n, loaded=0.0, num=30, nexp=1000000, cutoff=0.05, 
                    xmax=20, bins=50, figx=10, figy=8):
    '''
    n:      number of faces of die (or number of non-overlapping 
            events) (int)
    loaded: if 0.0, a fair die is presumed, otherwise loaded is a 
            float > 0: the larger the value, the more departure from a 
            uniform die distribution (float): it loaded is too large, 
            (it should be just a few tenths) negative p values may 
            result causing an error message
    '''
    die = 'loaded' if loaded > 0 else 'fair'
    probs = get_pvals(n, random=loaded)

    # the null hypothesis is the fair die, 
    # giving the expected number of each face
    null = num * get_pvals(n, random=0)

    # create observed values
    out = np.random.multinomial(num, probs, size=nexp)
    # form chi-square statistic
    stat = np.sum((out - null)**2 / null, axis=1)
    
    df   = n - 1 # degrees of freedom for chi-square
    tabl = create_chi_square_table()
    chi_value = tabl[df][1 - cutoff]
    prob_at_chi_value = chi2.pdf(chi_value, df=df)
    
    fig, (ax, ax1, ax2) = \
    plt.subplots(3, 1, gridspec_kw={'height_ratios': [5,1,1]}, 
                 figsize=(figx, figy))
    
    aa, bb, __ = ax.hist(stat, bins=bins, density=True)
    x = np.linspace(0, xmax, 1000)
    ax.plot(x, chi2.pdf(x, df), 'r', linewidth=4, 
            label='chi-square dist for DF: {}'.format(df))
    ax.set_xlim(0, xmax)
    ax.vlines(chi_value, color='k', linewidth=4, 
              ymin=0, ymax=prob_at_chi_value)
    ax.legend(fontsize='xx-large')
    txt = '{} of area is to the right\nof the black line'.format(cutoff)
    ax.text(chi_value, 1.5*prob_at_chi_value, txt, fontsize='xx-large')
    ax.grid()
    
    # null hypothesis p's histogram
    x = range(1, n + 1)
    ax1.bar(x, get_pvals(n, random=0), color='g', width=0.8)
    ax1.title.set_text('null hypothesis p values')
    ax1.grid()
    
    # data p's histogram
    ax2.bar(x, probs, color='c', width=0.8)
    ax2.title.set_text('data p values')
    ax2.grid
    
    txt =  'chi-square statistic for {}-sided {} die, '.format(n, die)
    txt += 'null hypothesis is fair die;'
    txt += '\nnumber of throws per trial: {}, number of trials: {}'\
           .format(num, nexp)
    fig.suptitle(txt, fontsize='xx-large')
    fig.tight_layout()
    plt.show()   

import pickle, gzip
# import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
# import seaborn as sns


def load_pgz(fname):
    with gzip.open(fname, 'rb') as f:
        o = pickle.load(f)
    return o

def _check_dim_error(v, ndim):
    if len(v.shape) != ndim:
        raise ValueError(f'require {ndim}-d array.')

def get_rank(vec):
    '''
    1-d array
    '''
    _check_dim_error(vec, 1)
    seq = np.arange(vec.shape[0])
    tmp = vec.argsort(axis=0) 
    rank = np.zeros((vec.shape[0]))
    rank[tmp] = seq
    return rank

def pval_obs2expected(pval_list):
    '''
    pval_list is list of m length-n vectors
    calculate column by column
    '''
    for i in pval_list:
        _check_dim_error(i, 1)
    o = []
    for i in pval_list:
        o.append((get_rank(i) + 1) / i.shape[0])
    return o

def qqplot(pval, ax, col_labels=None):
    if col_labels is not None and len(pval) != len(col_labels):
        raise ValueError('pval and col_labels has different number of instances.')
        
    pexp = pval_obs2expected(pval)
    n = len(pval)
    kw = {}
    for i in range(n):
        sort_idx = np.argsort(pexp[i])
        if col_labels is not None:
            kw = {'label': col_labels[i]}
        ax.scatter(-np.log(pexp[i][sort_idx]), -np.log(pval[i][sort_idx]), **kw)
    if col_labels is not None:
        ax.legend()

def p2chisq(p):
    return stats.norm.ppf(p / 2) ** 2

def get_median(gwas_code, gwas_dict, gwas_meta):
    gwas_codes = gwas_meta.gwas_id[gwas_meta.gwas_code == gwas_code].unique().tolist()
    nums = []
    for code in gwas_codes:
        nums += gwas_dict[code]['pvalue'].tolist()
    nums = np.array(nums)
    return np.median(p2chisq(nums))
    
def index_in_str_set(s, slist):
    for i in range(len(slist)):
        if s in slist[i]:
            return i, slist[i]
    return None, None

def z2p(z):
    return stats.norm.cdf(-np.abs(z)) * 2
    
import h5py
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import StrVector
pandas2ri.activate()
rbgen = importr('rbgen')

def query(pos, bgen_path, bgi_path, by='position'):
    if by == 'position':
        pos = [ int(p) for p in pos ]
        pos.sort()
        query = pd.DataFrame({
            'chromosome': [''],
            'start': [pos[0]], # int(gwas_df.pos[i])],
            'end': [pos[-1]], # [int(gwas_df.pos[i])]
        })
        kwarg = { 'ranges': query }
    elif by == 'id':
        query = StrVector(pos)
        kwarg = { 'rsids': query }
    else:
        raise ValueError('Args by is wrongly set.')

    print('querying', query)
    cached_data = rbgen.bgen_load(
        bgen_path,
        index_filename=bgi_path, 
        max_entries_per_sample=4,
        **kwarg
    )
    all_variants = pandas2ri.ri2py(cached_data[0])
    if all_variants.shape[0] != len(pos):
        raise ValueError('Extract fewer or more than input number of variants. Cannot handle.')
    all_probs = pandas2ri.ri2py(cached_data[4])
    return all_variants, all_probs

def query_hdf5(snp_start, nsnp, hdf5_path):
    with h5py.File(hdf5_path, 'r') as f:
        if f['position'].shape[0] <= snp_start + nsnp:
            raise ValueError('Too many SNPs: snp_start + nsnp >= # SNPs in the file.')
        geno = f['genotype'][:, snp_start : (nsnp + snp_start), :]
        pos = f['position'][snp_start : (nsnp + snp_start)]
    return pos, geno

def geno_prob_to_haplo(geno_prob):
    h1 = geno_prob[:, :, 1]
    h2 = geno_prob[:, :, 3]
    return h1, h2

def check_equal(mat1, mat2, thres=1e-6):
    if mat1.shape != mat2.shape:
        raise ValueError('Wrong dimension.')
    return np.absolute(mat1 - mat2) / mat1.shape[0] / mat1.shape[1] < thres


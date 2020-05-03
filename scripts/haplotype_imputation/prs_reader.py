import h5py
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '../prs')
import gwas_reader

def _get_samples(h5path):
    with h5py.File(h5path, 'r') as f:
        s = f['samples'][:].astype(str)
    return s

def load_prs(h5path, yaml):
    load_dict = gwas_reader.read_yaml(yaml)
    traits = list(load_dict.keys())
    indiv_list = _get_samples(h5path)
    prs = np.zeros((2, indiv_list.shape[0], len(traits)))
    # prs2 = np.array((len(indiv_list), len(traits)))
    
    with h5py.File(h5path, 'r') as f:
        trait_list = f['traits'][:].astype(str)
        pval_list = f['pval_cutoffs'][:]
        for idx, trait in enumerate(traits):
            trait_id = load_dict[trait]['name']
            pval = load_dict[trait]['pval']
            trait_idx = np.where(trait_list == trait_id)[0]
            pval_idx = np.where(pval_list == pval)[0]
            prs_mat = f['prs'][:, trait_idx, :, :][:, 0, pval_idx, :]
            for i in range(2):
                # breakpoint()
                prs[i, :, idx] = prs_mat[i, :]
    df_h = []
    for i in range(2):
        tmp_ = pd.DataFrame(prs[i, :, :])
        tmp_.columns = traits
        tmp_['individual_id'] = indiv_list
        df_h.append(tmp_)
    
    return df_h[0], df_h[1]
            
            
    


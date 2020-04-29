import argparse
parser = argparse.ArgumentParser(prog='run_haplo_logistic_solver.py', description='''
    Run logistic regression on a set of variants 
    in the mode of "haplo".
    It requires phenotype observes on both father and mother.
''')

parser.add_argument('--genotype-in-hdf5', help='''
    Genotype matrix in HDF5 format which is
    converted from BGEN by ../bgen2hdf5/ .
''')
parser.add_argument('--variant-chunk-size', type=int, default=50, help='''
    The number of variants to load & process together. 
''')
parser.add_argument('--n-threads', type=int, default=0, help='''
    The number of threads for reading HDF5 
    (default: no multi-threading at all).
''')
parser.add_argument('--gpu-index', type=int, default=0, help='''
    The index of GPU available
    (default = 0).
''')
parser.add_argument('--maf-filter', type=float, default=0.01, help='''
    Skip SNPs with MAF < maf-filter
    (default = 0.01).
    Range [0, 1].
''')
parser.add_argument('--father-phenotype-yaml', help='''
    YAML file specifying the phenotype of father.
    The structure follows:
        path: 'path_to_file'
        param:   # kwargs in pandas.read_csv
            header: 0
            sep: ','
        col:  # column name of phenotype of interest
            - 'alzheimer_disease'
            - 'heart_disease'
        indiv_col: 'eid'  # column name of individual id
''')
parser.add_argument('--mother-phenotype-yaml', help='''
    YAML file specifying the phenotype of mother.
    Follow the same structure as father phenotype YAML.
''')
parser.add_argument('--father-covariate-yaml', help='''
    ----------------------------------------------------
    [CAUTION]: Don't add sex. The software will assign 
    father as sex=0 and mother as sex=1 internally!
    ----------------------------------------------------
    YAML file specifying the covariates of father.
    Follow the same structure as father phenotype YAML.
''')
parser.add_argument('--mother-covariate-yaml', help='''
    ----------------------------------------------------
    [CAUTION]: Don't add sex. The software will assign 
    father as sex=0 and mother as sex=1 internally!
    ----------------------------------------------------
    YAML file specifying the covariates of mother.
    Follow the same structure as father phenotype YAML.
''')
parser.add_argument('--haplotype-imputation-yaml', default=None, help='''
    YAML file specifying the imputed haplotype
    parameterized as Pr(h1 is from father).
    Follow the same structure as father phenotype YAML.
''')
parser.add_argument('--out-npy', help='''
    Output in Python numpy NPY.
''')
args = parser.parse_args()

def _subset_row_by_indiv(df, list_):
    return df[ df['individual_id'].isin(list_) ].reset_index(drop=True)
def _extend_id(eid, par):
    return f'{eid}_{par}'
def _load_parent_in_pair(ff, mm, add_sex=False):
    # load the table
    pheno_df_dict = {
        'father': table_reader.load_table_from_yaml(ff),
        'mother': table_reader.load_table_from_yaml(mm),
    }
    # add sex an extra column encoding sex (M=0, F=1) 
    if add_sex is True:
        if 'sex' in pheno_df_dict['father'].columns:
            raise ValueError('sex should not occur in the table.')
        pheno_df_dict['father']['sex'] = 0
        pheno_df_dict['mother']['sex'] = 1
    # update the individual ID as individual ID + parent 
    return _add_parent_info_and_combine(pheno_df_dict)
def _extend_individual_id_list(id_list):
    '''
    Return pd data.frame with the extended individual ID and index. 
    '''
    dict_ = {}
    dict_['father'] = pd.DataFrame({
        'individual_id': id_list,
        'row_idx': [ i for i in range(len(id_list)) ]
    })   
    dict_['mother'] = dict_['father'].copy()
    return _add_parent_info_and_combine(dict_)
def _add_parent_info_and_combine(dict_):
    for parent in dict_.keys():
        dict_[parent]['individual_id'] = dict_[parent]['individual_id'].map(
            lambda x: _extend_id(x, parent)
        )
    df = pd.concat(
        (dict_['father'], dict_['mother']),
        axis=0
    )
    return df
def _load_haplotype_imputation(f):
    d_dict = {}
    d_dict['father'] = table_reader.load_table_from_yaml(f)
    d_dict['mother'] = d_dict['father'].copy()
    # for mother, Pr(h1 is from mother) = 1 - Pr(h1 is from father)
    d_dict['mother'].loc[:, d_dict['mother'].columns != 'individual_id'] = d_dict['mother'].loc[:, d_dict['mother'].columns != 'individual_id'].apply(lambda x: 1 - x)
    return _add_parent_info_and_combine(d_dict)
def multiply_mat_vec_col_by_col(mat, vec):
    return torch.einsum('ij,i->ij', mat, vec)

import logging, sys
# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

import time
import torch
from torch.utils import data
import pandas as pd
import numpy as np
from tqdm import tqdm
import table_reader
import geno_hdf5_reader
import data_stream
import logistic_gpu

if args.maf_filter > 1 or args.maf_filter < 0:
    raise ValueError(
        'maf should between 0 and 1, {} is given.'.format(args.maf_filter)
    )

logging.info('Loading phenotype')
pheno_df = _load_parent_in_pair(
    args.father_phenotype_yaml,
    args.mother_phenotype_yaml
)
logging.info('{} individual parents in phenotype'.format(pheno_df.shape[0]))

logging.info('Loading covariate')
covar_df = _load_parent_in_pair(
    args.father_covariate_yaml,
    args.mother_covariate_yaml,
    add_sex=True
)
logging.info('{} individual parents in covariate'.format(covar_df.shape[0]))

logging.info('Loading individual IDs from genotype')
individuals_in_genotype = geno_hdf5_reader.get_individual_list(args.genotype_in_hdf5)
# extend individual ID to individual ID + parent 
reference_indiv_df = _extend_individual_id_list(individuals_in_genotype)
logging.info('{} individual parents in genotype'.format(reference_indiv_df.shape[0]))

logging.info('Loading haplotype imputation result')
prob_z_df = _load_haplotype_imputation(args.haplotype_imputation_yaml)
logging.info('{} individual parents in haplotype imputation'.format(
    prob_z_df.shape[0])
)

logging.info('Jointing on individual parents')
# use individual id list from genotype as reference
# inner join all
reference_indiv_df = _subset_row_by_indiv(reference_indiv_df, pheno_df['individual_id'])
reference_indiv_df = _subset_row_by_indiv(reference_indiv_df, covar_df['individual_id'])
reference_indiv_df = _subset_row_by_indiv(reference_indiv_df, prob_z_df['individual_id'])

# Not yet implemented
# if args.individual_list is not None:
#     logging.info('-> Loading individual list')
#     indiv_list = table_reader.load_file_as_lines(args.individual_list)
#     indiv_extended_df = _extend_individual_id_list(indiv_list) 
#     logging.info('-> {} individuals in individual_list'.format(indiv_extended_df.shape[0]))
#     reference_indiv_df = _subset_row_by_indiv(reference_indiv_df, indiv_extended_df['individual_id'])

num_individual = reference_indiv_df.shape[0]
logging.info('{} individuals are extracted'.format(num_individual))

logging.info('Rearrange phenotype and covariate tables')
reference_pheno_mat = pd.merge(
    reference_indiv_df[['individual_id']], pheno_df, 
    left_on='individual_id', right_on='individual_id', how='left'
).drop(['individual_id'], axis=1)
reference_covar_mat = pd.merge(
    reference_indiv_df[['individual_id']], covar_df, 
    left_on='individual_id', right_on='individual_id', how='left'
).drop(['individual_id'], axis=1)
reference_prob_z_mat = pd.merge(
    reference_indiv_df[['individual_id']], prob_z_df, 
    left_on='individual_id', right_on='individual_id', how='left'
).drop(['individual_id'], axis=1)

num_phenotype = reference_pheno_mat.shape[1]
logging.info('{} phenotypes are extracted'.format(num_phenotype))
num_covariate = reference_covar_mat.shape[1]
logging.info('{} covariates are extracted'.format(num_covariate))
num_prob_z = reference_prob_z_mat.shape[1]
logging.info('{} prob z\'s are extracted'.format(num_prob_z))



logging.info('Prepare tensors')
column_y = reference_pheno_mat.columns.tolist()
y = torch.Tensor(reference_pheno_mat.values)
C = torch.cat(
    (
        torch.ones((num_individual, 1)), 
        torch.Tensor(reference_covar_mat.values)
    ),
    axis=1
)
z = torch.Tensor(reference_prob_z_mat.values)


logging.info('Prepare genotype stream')
params = {
    'batch_size': 1,
    'shuffle': False,
    'num_workers': args.n_threads
}
hdf5_reader = data_stream.HDF5ukbHapDataset(
    args.genotype_in_hdf5, 
    chunk_size=args.variant_chunk_size
)
variant_generator = data.DataLoader(hdf5_reader, **params)
num_variant = hdf5_reader.nvariant
logging.info('There are {} variants to work with'.format(num_variant))


logging.info('Run GWAS')
logging.info('-> Init solver')
solver = logistic_gpu.BatchLogisticSolver()
solver.message()

use_cuda = torch.cuda.is_available()
device = torch.device("cuda:{}".format(args.gpu_index) if use_cuda else "cpu")
y = y.to(device)
C = C.to(device)
z = z.to(device)

# out tensor: PrZ x Pheno x Variant x SumStat
out_tensor = torch.zeros((num_prob_z, num_phenotype, num_variant, 3)).to(device)

indiv_index = torch.LongTensor(reference_indiv_df['row_idx'].values)

niter = 0
snp_counter = 0
for h1, h2 in tqdm(variant_generator, total=hdf5_reader.nchunk):
    
    niter += 1
    
    # calculate step size and maf filter
    step_size = h1.shape[1]
    maf_filter = (
        h1.sum(axis=[0, 2]) + h2.sum(axis=[0, 2])
    ) / step_size / 2 > args.maf_filter
    
    # place holder for output
    bhat_ = torch.zeros((step_size,)).to(device)
    bse_ = torch.zeros((step_size, )).to(device)
    conv_ = torch.zeros((step_size, ), dtype=torch.bool).to(device)
    
    # haplotypes
    h1 = h1[0, :, indiv_index].T.to(device)
    h2 = h2[0, :, indiv_index].T.to(device)
    
    x_list = []
    for zi in range(num_prob_z):
        X = multiply_mat_vec_col_by_col(
            h1, z[:, zi]
        ) + multiply_mat_vec_col_by_col(
            h2, (1 - z[:, zi])
        )
        # X = X[indiv_index, :]
        X = X[:, maf_filter]
        x_list.append(X)
    X = torch.cat(x_list, axis=1)
    # t0 = time.time()
    # bhat, bse, conv = solver.batchIRLS(X.to(device), y[:, 0], C, device=device, use_mask=True, min_prob=1e-20)
    # t1 = time.time(); print('grand solve takes', t1 - t0)
    for p in range(num_phenotype):
        # t0 = time.time()
        bhat, bse, conv = solver.batchIRLS(X.to(device), y[:, p], C, device=device, use_mask=True, min_prob=1e-20)
        # t1 = time.time()
        bhat_[maf_filter] = bhat[-1]
        bse_[maf_filter] = bse[-1]
        conv_[maf_filter] = conv[-1]
        for zi in range(num_prob_z):
            start = zi * step_size
            end = (zi + 1) * step_size
            out_tensor[zi, p, snp_counter:(snp_counter + step_size), 0] = bhat_[start:end]
            out_tensor[zi, p, snp_counter:(snp_counter + step_size), 1] = bse_[start:end]
            out_tensor[zi, p, snp_counter:(snp_counter + step_size), 2] = conv_[start:end]
        # t2 = time.time(); print('solve takes', t1 - t0, ' save takes', t2 - t1)
    snp_counter += step_size
    if niter >= hdf5_reader.nchunk:
        break

out_tensor = out_tensor.to('cpu').numpy()
np.save(args.out_npy, out_tensor)


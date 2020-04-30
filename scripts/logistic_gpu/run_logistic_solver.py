import argparse
parser = argparse.ArgumentParser(prog='run_logistic_solver.py', description='''
    Run logistic regression on a set of variants.
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
parser.add_argument('--phenotype-yaml', help='''
    YAML file specifying the phenotype.
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
parser.add_argument('--covariate-yaml', help='''
    YAML file specifying the covariates.
    Follow the same structure as phenotype YAML.
''')
parser.add_argument('--individual-list', default=None, help='''
    The list of individual to be included in the GWAS
    (default = None which means to include all individuals
    who occur in both phenotype and covariate tables).
''')
parser.add_argument('--out-npy', help='''
    Output in Python numpy NPY.
''')
args = parser.parse_args()

def _subset_row_by_indiv(df, list_):
    return df[ df['individual_id'].isin(list_) ].reset_index(drop=True)

import logging, sys
# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

import torch
from torch.utils import data
import pandas as pd
import numpy as np
from tqdm import tqdm
import table_reader
import geno_hdf5_reader
import data_stream
import logistic_gpu

logging.info('Loading phenotype')
pheno_df = table_reader.load_table_from_yaml(args.phenotype_yaml)
logging.info('{} individuals in phenotype'.format(pheno_df.shape[0]))

logging.info('Loading covariate')
covar_df = table_reader.load_table_from_yaml(args.covariate_yaml)
logging.info('{} individuals in covariate'.format(covar_df.shape[0]))

logging.info('Loading individual IDs from genotype')
individuals_in_genotype = geno_hdf5_reader.get_individual_list(args.genotype_in_hdf5)
logging.info('{} individuals in genotype'.format(individuals_in_genotype.shape[0]))

logging.info('Jointing on individuals')
# use individual id list from genotype as reference
reference_indiv_df = pd.DataFrame({
    'individual_id': individuals_in_genotype,
    'row_idx': [ i for i in range(len(individuals_in_genotype)) ]
})     
# inner join all
reference_indiv_df = _subset_row_by_indiv(reference_indiv_df, pheno_df['individual_id'])
reference_indiv_df = _subset_row_by_indiv(reference_indiv_df, covar_df['individual_id'])

if args.individual_list is not None:
    logging.info('-> Loading individual list')
    indiv_list = table_reader.load_file_as_lines(args.individual_list)
    logging.info('-> {} individuals in individual_list'.format(len(indiv_list)))
    reference_indiv_df = _subset_row_by_indiv(reference_indiv_df, indiv_list)

num_individual = reference_indiv_df.shape[0]
logging.info('{} individuals are extracted'.format(num_individual))

logging.info('Rearrange phenotype and covariate tables')
reference_pheno_mat = pd.merge(reference_indiv_df[['individual_id']], pheno_df, left_on='individual_id', right_on='individual_id', how='left').drop(['individual_id'], axis=1)
reference_covar_mat = pd.merge(reference_indiv_df[['individual_id']], covar_df, left_on='individual_id', right_on='individual_id', how='left').drop(['individual_id'], axis=1)

num_phenotype = reference_pheno_mat.shape[1]
logging.info('{} phenotypes are extracted'.format(num_phenotype))
num_covariate = reference_covar_mat.shape[1]
logging.info('{} covariates are extracted'.format(num_covariate))


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
out_tensor = torch.zeros((3, num_phenotype, num_variant)).to(device)
indiv_index = torch.LongTensor(reference_indiv_df['row_idx'].values)
# np.save('cached_y.npy', y.to('cpu').numpy())
# np.save('cached_covar.npy', C.to('cpu').numpy())
# np.save('cached_indiv_index.npy', indiv_index.numpy())
niter = 0
snp_counter = 0
for h1, h2 in tqdm(variant_generator, total=hdf5_reader.nchunk):
    niter += 1
    X = h1[0, :, :].T + h2[0, :, :].T
    step_size = X.shape[1]
    maf_filter = X.sum(axis=0) / X.shape[0] / 2 > 0.01
    X = X[indiv_index, :]
    X = X[:, maf_filter]
    bhat_ = torch.zeros((step_size,)).to(device)
    bse_ = torch.zeros((step_size, )).to(device)
    conv_ = torch.zeros((step_size, ), dtype=torch.bool).to(device)
    for p in range(num_phenotype):
        bhat, bse, conv = solver.batchIRLS(X.to(device), y[:, p], C, device=device, use_mask=True, min_prob=1e-20)
        bhat_[maf_filter] = bhat[-1]
        bse_[maf_filter] = bse[-1]
        conv_[maf_filter] = conv[-1]
        out_tensor[0, p, snp_counter:(snp_counter + step_size)] = bhat_
        out_tensor[1, p, snp_counter:(snp_counter + step_size)] = bse_
        out_tensor[2, p, snp_counter:(snp_counter + step_size)] = conv_
    snp_counter += step_size
    if niter >= hdf5_reader.nchunk:
        break

out_tensor = out_tensor.to('cpu').numpy()
np.save(args.out_npy, out_tensor)


##
# Implement idea 2: multi-chromosome version
##

import argparse
parser = argparse.ArgumentParser(prog='impute_otf_multi_chr.py', description='''
    Impute parental origin of haplotypes and observed phenotypes.
    It takes preloaded phenotypes, covariates, and genotypes 
    (generated by impute_otf_preload.py)
''')

parser.add_argument('--genotype-prefix-pattern', help='''
    Prefix of preloaded NPY for genotype and position matrix. Should contain {chr_num} 
    as placeholder for chromosome number.
''')
parser.add_argument('--chromosomes', default=None, help='''
    List of chromosomes to work with, separated by ,.
    For instance: 1,2,3.
    If not set, it will include 1 .. 22.
''')
parser.add_argument('--npy-prefix', type=str, help='''
    Prefix of preloaded NPY for phenotypes and covariates.
''')
parser.add_argument('--output', help='''
    Output in TSV.GZ format.
''')
parser.add_argument('--nthread', default=None, type=int, help='''
    Number of threads to use.
''')
parser.add_argument('--imputer-output', type=str, help='''
    Pickle GZ imputer output
''')

args = parser.parse_args()


import logging, sys
import torch
import numpy as np
sys.path.insert(0, '../logistic_gpu')
import table_reader
import geno_hdf5_reader
import snp_list_reader
import haplotype_imputer
import gzip, pickle

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

if args.nthread is not None:
    torch.set_num_threads(args.nthread)

logging.info('Loading preloaded phenotypes and covariates')
    fmat = np.load(args.npy_prefix + '.fmat.npy')
    mmat = np.load(args.npy_prefix + '.mmat.npy')
    cmat = np.load(args.npy_prefix + '.cmat.npy')

logging.info('Loading posmat')
# load all posmat into memory since it does not take too much
if args.chromosomes is None:
    chroms = [ str(i) for i in range(1, 23) ]
else:
    chroms = args.chromosomes.split(',')
posmat_dic = {}
for chrom in chroms:
    posmat_dic[chrom] = np.load(args.genotype_prefix_pattern.format(chr_num=chrom) + '.posmat.npy')


logging.info('Run imputation: mode = multi-chromosome OTF'.format(args.impute_mode))
imputer = haplotype_imputer.HaploImputer()
beta, sigma2, out, lld = imputer.impute_otf_multi_chr(
    fmat, mmat, 
    posmat_dic, chroms, 
    args.genotype_prefix_pattern + '.hh1.npy', args.genotype_prefix_pattern + '.hh2.npy',
    cmat=cmat
)

logging.info('Save imputer output')
with gzip.open(args.imputer_output, 'w') as f:
    pickle.dump((beta, sigma2, lld), f)


logging.info('Output')
out.to_csv(args.output, compression='gzip', sep='\t', index=False)


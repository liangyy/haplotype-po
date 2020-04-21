import argparse
parser = argparse.ArgumentParser(prog='naive_prs.py', description='''
    Given:
        1. A list of GWAS summary statistic table;
        2. SNP map table (map each ukb_hap variant to GWAS variant);
        3. UKB HAP BGEN
        4. A list of P-value thresholds.
    Output:
        A N x K x P tensor where:
            1. N = sample size;
            2. K = number of GWASs;
            3. P = number of p-values.
''')
parser.add_argument('--gwas-yaml', help='''
    GWAS YAML:
        Each entry is a GWAS in the following struct:
        {
            gwas_name:
            {
                gwas:
                    {
                    path: path-to-table,
                    param:
                        {
                        compression: 'gzip'
                        }
                    }
                varcol: column name of variant ID (match with SNP map),
                betacol: column name of effect size,
                pvalcol: column name of p-value,
                ld_clump: path to LD clump table (optional)
            }
        }
''')
parser.add_argument('--snp-map', help='''
    SNP map
''')
parser.add_argument('--bgen', help='''
    BGEN file path
''')
parser.add_argument('--bgi', help='''
    BGEN BGI file path
''')
parser.add_argument('--sample', help='''
    BGEN SAMPLE file path
''')
parser.add_argument('--pval-cutoffs', help='''
    P-value cutoffs in naive PRS.
    Separated by ','
''')
parser.add_argument('--chromosome', type=str, help='''
    The chromosome number we are working with
    (account for the missing chromosome in ukb_hap_v2)
''')
parser.add_argument('--output-hdf5', help='''
    HDF5 file name of output (if not exists, it will be created)
''')
parser.add_argument('--snp-chunk-size', type=int, default=100, help='''
    Number of SNPs to process at a time
''')
parser.add_argument('--bgen-writing-cache-size', type=int, default=50, help='''
    BGEN reading cache size in MB.
''')
parser.add_argument('--max-sample-chunk-size', type=int, default=10000, help='''
    Maximum size of chunk on sample axis.
''')
parser.add_argument('--max-trait-chunk-size', type=int, default=5, help='''
    Maximum size of chunk on trait axis. 
''')
args = parser.parse_args()

import logging, time, sys
# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

from tqdm import tqdm
import pandas as pd
import gwas_reader
import ukb_hap_reader
import prs_matrix

def convert_pval_args_str(instr):
    try:
        return [ float(i) for i in instr.split(',') ]
    except:
        raise ValueError('Wrong pval-cutoffs: {instr}')

def size_in_mb_to_cache_size(size):
    return int(size * (1024 ** 2))    

pval_cutoffs = convert_pval_args_str(args.pval_cutoffs)

logging.info('Loading GWAS')
gwas_dict = gwas_reader.gwas_reader(
    args.gwas_yaml, 
    snp_map=args.snp_map,
    logger=logging
)

logging.info('Generating variant list')
var_df = gwas_reader.build_var_df(gwas_dict, logger=logging)

logging.info('Build BGEN reader')
hap_reader = ukb_hap_reader.UKBhapReader(
    bgen_path=args.bgen, 
    bgen_bgi_path=args.bgi, 
    sample_path=args.sample
)
var_generator = hap_reader.retrieve_from_list(
    # NOTE
    # chrom is empty due to ukb_hap_v2 has missing chromosome
    chrom=[ '' for i in range(var_df.shape[0]) ], 
    pos=var_df.pos.tolist(), 
    non_effect_allele=var_df.effect_allele.tolist(), 
    effect_allele=var_df.non_effect_allele.tolist(),
    n_var_cached=args.snp_chunk_size
)

logging.info('Initialize PRS matrix')
prs_mat = prs_matrix.PRSmatrix(
    gwas_dict, args.sample, args.chromosome, pval_cutoffs, args.output_hdf5,
    cache_size=size_in_mb_to_cache_size(args.bgen_writing_cache_size), 
    max_sample_chunk_size=args.max_sample_chunk_size, 
    max_trait_chunk_size=args.max_trait_chunk_size
)

logging.info('Update PRS')
for dosage_row in tqdm(var_generator, total=var_df.shape[0]):
    prs_mat.update(dosage_row)

logging.info('Save PRS')
prs_mat.save(logger=logging)


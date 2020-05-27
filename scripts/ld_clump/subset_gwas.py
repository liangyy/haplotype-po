import argparse
parser = argparse.ArgumentParser(prog='subset_gwas.py', description='''
    Subset GWAS file so that it only contains SNPs in 
    SNP map
''')
parser.add_argument('--gwas', help='''
    Input GWAS file name (assume tsv.gz).
''')
parser.add_argument('--snp-map', help='''
    The SNP map generated in ../prepare_snp_map/.
    Assume tsv.gz.
''')
parser.add_argument('--output', help='''
    File name of output.
''')
parser.add_argument('--gwas-variant-col', help='''
    column name of variant ID.
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

import pandas as pd

logging.info('Loading GWAS')
gwas = pd.read_csv(
    args.gwas, 
    header=0, 
    sep='\t', 
    compression='gzip'
)

logging.info('Loading SNP map')
snp_map = pd.read_csv(
    args.snp_map,
    sep='\t',
    header=0,
    compression='gzip',
    dtype={3:'str'}
)

logging.info('Intersect')
snp_map = snp_map[ (snp_map['assigned_id'] != 'not_shown') & (snp_map['assigned_id'] != 'ambiguious') ]
gwas_sub = gwas[ gwas[args.gwas_variant_col].isin(snp_map['assigned_id']) ]
gwas_cols = gwas_sub.columns.tolist()
gwas_sub = pd.merge(gwas_sub, snp_map[['chrom', 'pos', 'allele_ids', 'assigned_id']], left_on='variant', right_on='assigned_id')
gwas_sub['snp_id_in_ukb_hap'] = gwas_sub['assigned_id']

# save 
logging.info('Saving the results')
gwas_sub = gwas_sub[ gwas_cols + [ 'snp_id_in_ukb_hap' ] ]
gwas_sub.to_csv(args.output, compression='gzip', sep='\t', index = None)


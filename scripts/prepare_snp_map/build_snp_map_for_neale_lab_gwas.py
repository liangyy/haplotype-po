import argparse
parser = argparse.ArgumentParser(prog='build_snp_map_for_neale_lab_gwas.py', description='''
    Build the SNP map table: phased genotype variant <=> Neale's lab GWAS
''')
parser.add_argument('--genotype-pattern', help='''
    In the form: prefix{chr}suffix.
    Will load 1..22 chromosomes (no X).
''')
parser.add_argument('--genotype-sample', help='''
    The corresponding sample file
''')
parser.add_argument('--output', help='''
    File name of output (if not exists, it will be created)
''')
parser.add_argument('--gwas', help='''
    Neale's lab GWAS (one GWAS as example, 
    they all share the same variant set)
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

import bgen_reader
import pandas as pd
import helper

logging.info('Loading GWAS')
gwas = pd.read_csv(args.gwas, header=0, sep= '\t', compression='gzip')
map_table = pd.DataFrame()
for i in range(1, 23):
    i = str(i)
    
    logging.info(f'Processing chr{i}: Loading BGEN')
    bgen = bgen_reader.read_bgen(
        args.genotype_pattern.format(chr=i), 
        samples_filepath = args.genotype_sample
    )
    
    logging.info(f'Processing chr{i}: Loading variant table')
    variant = bgen["variants"].compute()
    variant['chrom'] = i
    
    logging.info(f'Processing chr{i}: Building variant ID candidates')
    variant['allele_1st'] = variant['allele_ids'].apply(lambda x: x.split(',')[0])
    variant['allele_2nd'] = variant['allele_ids'].apply(lambda x: x.split(',')[1])
    variant['varid1'] = variant[['chrom', 'pos', 'allele_1st', 'allele_2nd']].apply(lambda x: helper.make_id(x), axis=1)
    variant['varid2'] = variant[['chrom', 'pos', 'allele_2nd', 'allele_1st']].apply(lambda x: helper.make_id(x), axis=1)
    
    logging.info(f'Processing chr{i}: Running checker')
    variant_check = helper.join_with_varid(
        variant['varid1'], 
        variant['varid2'], 
        gwas['variant']
    )
    variant = pd.merge(variant, variant_check, left_on=['varid1', 'varid2'], right_on=['id1', 'id2'], how='left')
    
    map_table = pd.concat([map_table, variant[['chrom', 'pos', 'allele_ids', 'id', 'rsid', 'assigned_id']]])
    
# save 
logging.info('Saving the results')
map_table.to_csv(args.output, compression='gzip', sep='\t', index = None)


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
parser.add_argument('--GWAS-yaml', help='''
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
    P-value cutoffs in naive PRS
''')
parser.add_argument('--output-hdf5', help='''
    HDF5 file name of output (if not exists, it will be created)
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
    
    map_table = pd.concat([map_table, variant[['chrom', 'pos', 'allele_ids', 'id', 'rsid', 'assigned_id', 'assigned_sign']]])
    
# save 
logging.info('Saving the results')
map_table.to_csv(args.output, compression='gzip', sep='\t', index = None)


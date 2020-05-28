import argparse
parser = argparse.ArgumentParser(prog='annot_rsid.py', description='''
    Add rsID to Neale lab UKBB GWAS
''')

parser.add_argument('--gwas', help='''
    Neale lab UKBB GWAS summary statistics
''')
parser.add_argument('--variant-meta', help='''
    Meta-information of the variants in Neale lab UKBB GWAS
''')
parser.add_argument('--output', help='''
    Output
''')
args = parser.parse_args()
import pandas as pd

gwas = pd.read_csv(args.gwas, sep='\t', compression='gzip')
var_meta = pd.read_csv(args.variant_meta, sep='\t', compression='gzip')
gwas = pd.merge(gwas, var_meta[['variant', 'rsid', 'ref', 'alt']], left_on='variant', right_on='variant', how='inner')
gwas.rename(
    columns = {
        'ref':'non_effect_allele',
        'alt': 'effect_allele'
    }, 
    inplace = True
)

gwas.to_csv(args.output, sep='\t', compression='gzip', index=False)


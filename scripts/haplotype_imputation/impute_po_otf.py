##
# Implement idea 2
##

import argparse
parser = argparse.ArgumentParser(prog='impute_po_otf.py', description='''
    Impute parent of origin from haplotypes and observed phenotypes.
    For each trait, provide a list of variants to work with.
''')

parser.add_argument('--genotype-in-hdf5', help='''
    Genotype matrix in HDF5 format which is
    converted from BGEN by ../bgen2hdf5/ .
''')
parser.add_argument('--snp-list-yaml', help='''
    YAML to specify the list of SNPs to work with. 
    Example:
    trait1:
        load_method: 'ld_clump'
        params:
            gwas_yaml: 
                gwas: 
                    path: 'test_inputs/test_gwas1.txt'
                    param:
                      header: 0
                      sep: ' '
                varcol: 'variant'
                betacol: 'beta'
                pvalcol: 'pval'
                ld_clump: 'test_inputs/test_gwas1.clump'
            snp_map: 'path_to_snp_map'
            pval_cutoff: 0.05
    trait2:
        load_method: 'from_list'
        params:
            list_path: 'path_to_list'
''')
parser.add_argument('--snp-list-cache', type=str, help='''
    Path to PGZ file (pickle.gz) file to cache 
    the information loaded from --snp-list-yaml.
    Require pgz as extension.
''')
parser.add_argument('--father-phenotype-yaml', help='''
    Observed phenotype of father
    with all loading specified in YAML:
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
    Observed phenotype of mother 
    with all loading specified in YAML 
    (follow the same structure as father's)
''')
parser.add_argument('--shared-covariate-yaml', default=None, help='''
    Covariate shared by father and mother
    with all loading specified in YAML 
    (follow the same structure as phenotype's)
''')
parser.add_argument('--impute-mode', type=str, default='basic_em', help='''
    Imputation approach.
''')
parser.add_argument('--output', help='''
    Output in TSV.GZ format.
''')
args = parser.parse_args()


import logging, sys
sys.path.insert(0, '../logistic_gpu')
import table_reader
import geno_hdf5_reader
import snp_list_reader
import haplotype_imputer

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)



logging.info('Loading observed phenotypes')
df_father = table_reader.load_table_from_yaml(
    args.father_phenotype_yaml,
    rename_cols=True
)
df_mother = table_reader.load_table_from_yaml(
    args.mother_phenotype_yaml,
    rename_cols=True
)

df_covar = None
if args.shared_covariate_yaml is not None:
    logging.info('Loading covariates')
    df_covar = table_reader.load_table_from_yaml(
        args.shared_covariate_yaml,
        rename_cols=True
    )
    # breakpoint()
    df_covar = table_reader.standardize_columns(df_covar, except_cols=['individual_id'])

logging.info('Loading variant list')
snp_loader = snp_list_reader.snpLoader(args.snp_list_yaml)
snp_loader.load(args.snp_list_cache)

logging.info('Loading haplotypes')
h1, h2, hap_indiv_df, hap_pos_df = geno_hdf5_reader.load_haplotypes_by_position(
    args.genotype_in_hdf5, snp_loader.snp_pos_dict
)
# breakpoint()

logging.info('Run imputation: mode = {}'.format(args.impute_mode))
imputer = haplotype_imputer.HaploImputer()
_, _, out, lld = imputer.impute_otf(
    df_father, df_mother, 
    h1, h2, hap_indiv_df, hap_pos_df,
    mode=args.impute_mode,
    df_covar=df_covar
)
# lld_ = [ l.numpy()[0] for l in lld ]
# logging.info('lld = ', ' '.join(lld))
print(lld)

logging.info('Output')
out.to_csv(args.output, compression='gzip', sep='\t', index=False)


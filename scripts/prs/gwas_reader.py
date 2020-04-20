import yaml
import pandas as pd

def read_yaml(yaml):
    with open(yaml) as file:
        dict = yaml.load(file)
    return dict
def read_snp_map(path):
    if path is None:
        return None
    snp_map_df = pd.read_csv(snp_map)
    if 'allele_ids' in snp_map_df.columns:
        snp_map_df['effect_allele'] = snp_map_df['allele_ids'].map(lambda x: x.split(',')[0])
        snp_map_df['non_effect_allele'] = snp_map_df['allele_ids'].map(lambda x: x.split(',')[1])
    return snp_map_df
def gwas_reader(yaml, snp_map=None, logger=None):
    gwas_dict = read_yaml(yaml)
    snp_map_df = read_snp_map(snp_map)
    out_dict = {}
    counter = 0
    ntotal = len(gwas_dict.keys())
    for i in gwas_dict.keys():
        
        # print some log
        counter += 1
        message = f'gwas_reader: processing {i}, {counter}/{ntotal}'
        if logger is not None:
            logger.info(message)
        else:
            print(message)
        # END
        
        this_gwas = gwas_dict[i]
        gwas_df = pd.read_csv(
            this_gwas['gwas']['path'], 
            usecols=[
                this_gwas['varcol'],
                this_gwas['betacol'],
                this_gwas['pvalcol']
            ],
            names=[
                'variant_id',
                'effect_size',
                'pvalue'
            ]
            **this_gwas['gwas']['param']
        )
        
        if 'ld_clump' in this_gwas:
            clump_df = pd.read_csv(
                this_gwas['ld_clump'], 
                header=None,
                names=['variant_id']
            )
            gwas_df = gwas_df[ 
                gwas_df['variant_id'].isin(clump_df['variant_id']) 
            ]
        
        if snp_map_df is None:
            pass
        else:
            gwas_df = pd.merge(
                gwas_df, 
                snp_map_df,
                how='left',
                left_on=['variant_id'],
                right_on=['assigned_id'],
                suffixes=('_gwas', '_snp_map')
            )
        out_dict[i] = gwas_df
    return out_dict
            
        
        
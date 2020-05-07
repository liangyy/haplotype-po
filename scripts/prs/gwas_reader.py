import yaml, re, os
import pickle, gzip
import pandas as pd

def read_yaml(yamlfile):
    with open(yamlfile) as f:
        mydict = yaml.safe_load(f)
    return mydict

def _log_print(msg, logger):
    if logger is None:
        print(msg)
    else:
        logger.info(msg)

def read_snp_map(path):
    if path is None:
        return None
    snp_map_df = pd.read_csv(
        path, 
        compression='gzip', 
        sep='\t', 
        header=0,
        dtype={3:'str'}
    )
    if 'allele_ids' in snp_map_df.columns:
        snp_map_df['effect_allele'] = snp_map_df['allele_ids'].map(lambda x: x.split(',')[0])
        snp_map_df['non_effect_allele'] = snp_map_df['allele_ids'].map(lambda x: x.split(',')[1])
    return snp_map_df

def clean_tab(mydict):
    for i in mydict.keys():
        if isinstance(mydict[i], str) and '\\t' in mydict[i]:
            mydict[i] = re.sub(mydict[i], '\\t', '\t')
    return mydict

def gwas_reader(yaml, snp_map=None, logger=None, cache_path=None, from_dict=True):
    if cache_path is not None and os.path.exists(cache_path) and os.path.isfile(cache_path):
        filename, file_extension = os.path.splitext(cache_path)
        desired_ext = '.pgz'
        if file_extension != desired_ext:
            raise ValueError(f'cache_path should have {desired_ext} as extension.')
        _log_print(f'Loading gwas_dict cached in {cache_path}', logger)
        with gzip.open(cache_path, 'rb') as f:
            out_dict  = pickle.load(f)
        return out_dict
    if from_dict is False:
        gwas_dict = read_yaml(yaml)
    else:
        gwas_dict = yaml
    snp_map_df = read_snp_map(snp_map)
    out_dict = {}
    counter = 0
    ntotal = len(gwas_dict.keys())
    for i in gwas_dict.keys():
        
        # print some log
        counter += 1
        message = f'gwas_reader: processing {i}, {counter}/{ntotal}'
        _log_print(message, logger)
        # END
        
        this_gwas = gwas_dict[i]
        if 'param' not in this_gwas['gwas']:
            this_gwas['gwas']['param'] = {}
        else:
            this_gwas['gwas']['param'] = clean_tab(this_gwas['gwas']['param'])
        
        this_gwas['gwas']
        gwas_df = pd.read_csv(
            this_gwas['gwas']['path'], 
            usecols=[
                this_gwas['varcol'],
                this_gwas['betacol'],
                this_gwas['pvalcol']
            ],
            **this_gwas['gwas']['param']
        )
        gwas_df.columns = ['variant_id', 'effect_size', 'pvalue']
        
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
        
        # add my_var_id in the format: chr:pos:ref:alt
        gwas_df['my_var_id'] = gwas_df[['chrom', 'pos', 'effect_allele', 'non_effect_allele']].apply(
            lambda x: '{}:{}:{}:{}'.format(x.chrom, x.pos, x.effect_allele, x.non_effect_allele),
            axis=1
        )
        
        out_dict[i] = gwas_df
    if cache_path is not None:
        dirname = os.path.dirname(cache_path)
        if not os.path.exists(dirname):
            raise ValueError('Expect the directory to cache_path exists.')
        else:
            _log_print(f'Caching gwas_dict to {cache_path}', logger)
            with gzip.open(cache_path, 'wb') as f:
                pickle.dump(out_dict, f)
    return out_dict
            
def build_var_df(gwas_dict, logger=None):
    '''
    aggregate all variants in one data.frame
    '''    
    # lazy load, not really using the reader 
    var_df = pd.DataFrame()
    for i in gwas_dict.keys():
        
        # some log message
        message = f'build_var_df: processing {i}'
        _log_print(message, logger)
        # END
        
        tmp_ = gwas_dict[i][['chrom', 'pos', 'effect_allele', 'non_effect_allele']]
        var_df = pd.concat([var_df, tmp_], axis=0).drop_duplicates(keep='first').reset_index(drop=True)
        
    return var_df
        

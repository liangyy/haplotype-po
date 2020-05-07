import sys
sys.path.insert(0, '../prs')
import gwas_reader
import os
import pickle, gzip


class snpLoader:
    def __init__(self, snp_yaml):
        self.snp_yaml = snp_yaml
        self._load_yaml()
    @staticmethod
    def _check_element(dict_, ele_):
        if ele_ not in dict_:
            raise ValueError(f'Missing {ele_}.')
    def load(self, cache_path=None):
        
        # load cached file in there is any
        if cache_path is not None and os.path.exists(cache_path) and os.path.isfile(cache_path):
            filename, file_extension = os.path.splitext(cache_path)
            desired_ext = '.pgz'
            if file_extension != desired_ext:
                raise ValueError(f'cache_path should have {desired_ext} as extension.')
            print(f'Loading gwas_dict cached in {cache_path}')
            with gzip.open(cache_path, 'rb') as f:
                self.snp_pos_dict = pickle.load(f)
            return
        
        # load from scratch is no cached file
        self.snp_pos_dict = {}
        for trait in self.loader_dict.keys():
            load_method = getattr(
                self, 
                '_load_{}'.format(self.loader_dict[trait]['load_method'])
            )
            self.snp_pos_dict[trait] = load_method(
                self.loader_dict[trait]['params']
            )
        
        # cache self.loader_dict if having cache_path
        if cache_path is not None:
            dirname = os.path.dirname(cache_path)
            if not os.path.exists(dirname):
                raise ValueError('Expect the directory to cache_path exists.')
            else:
                print(f'Caching gwas_dict to {cache_path}')
                with gzip.open(cache_path, 'wb') as f:
                    pickle.dump(self.snp_pos_dict, f)
                     
    def _load_yaml(self):
        snp_yaml_dict = gwas_reader.read_yaml(self.snp_yaml)
        self.loader_dict = {}
        for trait in snp_yaml_dict.keys():
            self._check_element(snp_yaml_dict[trait], 'load_method')
            self._check_element(snp_yaml_dict[trait], 'params')
            self.loader_dict[trait] = snp_yaml_dict[trait]
    def _load_ld_clump(self, param_dict):
        self._check_element(param_dict, 'gwas_yaml')
        self._check_element(param_dict, 'snp_map')
        self._check_element(param_dict, 'pval_cutoff')
        gwas_yaml = {}
        gwas_yaml['gwas_yaml'] = param_dict['gwas_yaml']
        out = self.__load_ld_clump(
            gwas_yaml=gwas_yaml,
            snp_map=param_dict['snp_map'],
            pval_cutoff=param_dict['pval_cutoff']
        )
        return out
    def __load_ld_clump(self, gwas_yaml, snp_map, pval_cutoff):
        df_gwas = gwas_reader.gwas_reader(gwas_yaml, snp_map=snp_map, logger=None, cache_path=None, from_dict=True)['gwas_yaml']
        df_gwas = df_gwas[ df_gwas['pvalue'] < pval_cutoff ]
        pos = df_gwas['pos'].tolist()
        return pos
        
        

def load(snp_yaml):
    
        
    return 

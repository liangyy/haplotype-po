import sys
sys.path.insert(0, '../prs')
import gwas_reader

class snpLoader:
    def __init__(self, snp_yaml):
        self.snp_yaml
        self._load_yaml()
    @staticmethod
    def _check_element(dict_, ele_):
        if ele_ not in dict_:
            raise ValueError(f'Missing {ele_}.')
    def load(self):
        self.snp_pos_dict = {}
        for trait in self.loader_dict.keys():
            load_method = getattr(
                self, 
                f'_load_{}'.format(trait['load_method'])
            )
            self.snp_pos_dict[trait] = load_method(
                loader_dict['params']
            )
    def _load_yaml(self):
        snp_yaml_dict = gwas_reader.read_yaml(snp_yaml)
        self.loader_dict = {}
        for trait in snp_yaml_dict.keys():
            self._check_element(self.loader_dict, 'load_method')
            self._check_element(self.loader_dict, 'params')
            self.loader_dict[trait] = snp_yaml_dict[trait]
    def _load_ld_clump(self, param_dict):
        self._check_element(param_dict, 'gwas_yaml')
        self._check_element(param_dict, 'snp_map')
        self._check_element(param_dict, 'pval_cutoff')
        gwas_yaml = {}
        gwas_yaml['gwas_yaml'] = param_dict['gwas_yaml']
        out = __load_ld_clump(
            gwas_yaml=gwas_yaml,
            snp_map=param_dict['snp_map']
            pval_cutoff=param_dict['pval_cutoff']
        )
        return out
    def __load_ld_clump(self, gwas_yaml, snp_map, pval_cutoff):
        df_gwas = gwas_reader.gwas_reader(gwas_yaml, snp_map=snp_map, logger=logger, cache_path=None)['gwas_yaml']
        df_gwas = df_gwas[ df_gwas['pvalue'] < pval_cutoff ]
        pos = df_gwas['pos'].tolist()
        return pos
        
        

def load(snp_yaml):
    
        
    return 
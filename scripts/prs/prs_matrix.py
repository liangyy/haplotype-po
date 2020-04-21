import h5py
import numpy as np

# reference: https://github.com/liangyy/predixcan_prediction/blob/a85d52d89de9fe237a1217b5627c7e8d9f700f7e/predict.py#L70
class PRSmatrix:
    def __init__(self, gwas_dict, bgen_sample, chromosome, pval_cutoffs, output_h5, 
        cache_size=int(50 * (1024 ** 2)), 
        max_sample_chunk_size=10000, # this has been optimized for UKB sample size
        max_trait_chunk_size=5):
        self.H5_file = None
        self.gwas_dict = gwas_dict
        self.gwas_index = { k: i for i, k in enumerate(gwas_dict) }
        self.bgen_sample = bgen_sample
        self.chromosome = chromosome
        self.output_h5 = output_h5
        self.cache_size = cache_size
        self.max_sample_chunk_size = max_sample_chunk_size
        self.max_trait_chunk_size = max_trait_chunk_size
        self.nsample = self._get_n_sample()
        self.pval_cutoffs = np.array(pval_cutoffs)
        self.nhap = 2  # NOTE: assume diploid
        self._index_variant()
    
    def _index_variant(self):
        if hasattr(self, 'gwas_variant_index'):
            return
        else:
            self.gwas_variant_index = {}
            for i in self.gwas_index.keys():
                inner_dict = { self.gwas_dict[i].my_var_id[i_] : i_ for i_ in range(self.gwas_dict[i].shape[0]) }
                self.gwas_variant_index[i] = inner_dict
    
    @property
    def ntrait(self):
        return len(self.gwas_dict)
        
    @property
    def ncutoff(self):
        return self.pval_cutoffs.shape[0]
    
    def _var_in(self, gwas_name, my_var_id):
        if gwas_name not in self.gwas_variant_index:
            raise ValueError(f'{gwas_name} is not in gwas_variant_index/gwas_dict.')
        return my_var_id in self.gwas_variant_index[gwas_name]
    
    def _get_gwas_info_by_var_id(self, gwas_name, info_col, var_id):
        '''
        If var_id not in gwas_name, return None
        '''
        if not self._var_in(gwas_name, var_id):
            return None
        else:
            row_idx_of_var_id = self.gwas_variant_index[gwas_name]['var_id']
            return self.gwas_dict[gwas_name][info_col][row_idx_of_var_id]
        
    def _get_n_sample(self):
        with open(self.bgen_sample, 'r') as samples:
            counter = 0
            for line_idx, line in enumerate(samples):
                if line_idx <= 1:
                    continue
                counter += 1
        return counter
    
    def _get_samples(self):
        with open(self.bgen_sample, 'r') as samples:
            for line_idx, line in enumerate(samples):
                if line_idx <= 1:
                    continue
                line_split = line.split()
                yield [line_split[0], line_split[1]]

    @staticmethod
    def _get_effect_size(beta, sign):
        if sign == '+':
            return beta
        elif sign == '-':
            return -beta
        else:
            raise ValueError(f'sign = {sign} is not allowed.')
    
    def update(self, dosage_row):
        if self.H5_file is None:
            self.H5_file = h5py.File(
                self.output_h5, 'w', 
                # chunk_cache_mem_size=self.cache_size
                rdcc_nbytes=self.cache_size
            )
            
            if self.nsample != len(dosage_row.haplo_dosage_1):
                raise ValueError('the length of sample file ({}) is different from the BGEN file ({})'.format(
                    self.nsample,
                    len(dosage_row.haplo_dosage_1)
                ))
            
            n_trait_chunk = np.min((self.ntrait, self.max_trait_chunk_size))
            n_sample_chunk = np.min((self.nsample, self.max_sample_chunk_size))
            n_cutoff_chunk = self.ncutoff
            n_hap_chunk = self.nhap
            self.H5_prs = self.H5_file.create_dataset(
                "prs", 
                shape=(self.nsample, self.ntrait, self.ncutoff, self.nhap), # sample x trait x cutoff x hap
                chunks=(n_sample_chunk, n_trait_chunk, n_cutoff_chunk, n_hap_chunk),
                dtype=np.dtype('float32'), 
                scaleoffset=4, 
                compression='gzip'
            )
        # traverse all gwas results and all p-value cutoffs
        for i in self.gwas_index.keys():
            pval = self._get_gwas_info_by_var_id(i, 'pvalue', dosage_row.my_var_id)
            if pval is None:
                continue
            beta = self._get_gwas_info_by_var_id(i, 'effect_size', dosage_row.my_var_id)
            sign = self._get_gwas_info_by_var_id(i, 'assigned_sign', dosage_row.my_var_id)
            effect_size = self._get_effect_size(beta, sign)
            gwas_idx = self.gwas_index[i]
            self.H5_prs[:, gwas_idx, pval < self.pval_cutoffs, 0] = dosage_row.haplo_dosage_1 * effect_size
            self.H5_prs[:, gwas_idx, pval < self.pval_cutoffs, 1] = dosage_row.haplo_dosage_2 * effect_size
                # yield i, dosage_row
    
    @staticmethod
    def _myprint(msg, logger):
        if logger is None:
            print(msg)
        else:
            logger.info(msg)
    
    def save(self, logger=None):
        
        # samples
        sample_generator = self._get_samples()
        self.H5_sample = self.H5_file.create_dataset("samples", (self.nsample,), dtype='S25')
        for col in range(0, self.H5_prs.shape[1]):
            try:
                self.H5_sample[col] = np.string_(next(sample_generator)[0])
            except StopIteration:
                self._myprint("ERROR: There are not enough rows in your sample file!", logger)
                self._myprint(
                    "Make sure dosage files and sample files have the same number of individuals in the same order.", 
                    logger
                )
                os.remove(self.output_h5)
                sys.exit(1)
        
        # gwas's
        gwas_dset = self.H5_file.create_dataset("traits", (self.ntrait,), dtype='S30')
        for gwas_name in self.gwas_index.keys():
            gwas_dset[gwas_index[gwas_name]] = np.string_(str(gwas_name))
        
        # p-values
        pval_dset = self.H5_file.create_dataset("pval_cutoffs", (self.ncutoff,), dtype='float32')
        pval_dset = self.pval_cutoffs

        self.H5_file.close()

        # check number of samples
        try:
            next(sample_generator)
        except StopIteration:
            self._myprint(
                "{} Predicted expression file complete!".format(datetime.datetime.now()),
                logger
            )
        else:
            self._myprint("ERROR: There are too many rows in your sample file!", logger)
            self._myprint(
                "Make sure dosage files and sample files have the same number of individuals in the ame order.", 
                logger
            )
            if os.path.isfile(self.output_h5):
                os.remove(self.output_h5)
            sys.exit(1)        
                                                
                                                
                                                

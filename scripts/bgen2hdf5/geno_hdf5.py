import h5py
import numpy as np

class GenotypeHDF5:
    def __init__(list_sample, num_variants, output_h5, dict_variant_meta,
        phased=True, dtype=int, chromosome=None, cache_size=int(50 * (1024 ** 2),
        max_sample_chunk_size=10000, max_variant_chunk_size=100):
        self.H5_file = None
        self.ARRAY_genotype = None
        self.list_sample = list_sample
        self.nvariant = num_variants
        self.chromosome = chromosome
        self.output_h5 = output_h5
        self.cache_size = cache_size
        self.max_sample_chunk_size = max_sample_chunk_size
        self.max_variant_chunk_size = max_variant_chunk_size
        self.phased = phased
        self.dtype = dtype 
        self.n_filled_variants = 0 
        self.dict_variant_meta = self._init_dict_variant_meta(dict_variant_meta)
    
    def _init_dict_variant_meta(self, dict_variant_meta):
        meta_pool = ('position', 'chromosome', 'reference_allele', 'dosage_allele')
        self.dict_variant_meta = {}
        for i in dict_variant_meta.keys():
            if i in meta_pool:
                self.dict_variant_meta[i] = (
                    dict_variant_meta[i],
                    []
                )
    
    @property
    def nhap(self):
        if self.phased is True:
            return 2
        else:
            return 1
    
    @property
    def nsample(self):
        return len(self.list_sample)
    
    def fill_in(self, dosage_row, genotype_entries):
        
        if len(genotype_entries) != self.nhap:
            raise ValueError(f'genotype_entries = {genotype_entries} but we set nhap with a different number.')
        
        if self.H5_file is None:
            self.H5_file = h5py.File(
                self.output_h5, 'w', 
                rdcc_nbytes=self.cache_size
            )
            size_variant_chunk = np.min((self.nvariant, self.max_variant_chunk_size))
            size_sample_chunk = np.min((self.nsample, self.max_sample_chunk_size))
            size_hap_chunk = 1
            self.ARRAY_genotype = self.H5_file.create_dataset(
                "genotype", 
                # hap x variant x sample (C order)
                shape=(self.nhap, self.nvariant, self.nsample), 
                chunks=(size_hap_chunk, size_variant_chunk, size_sample_chunk),
                dtype=self.dtype, 
                compression='gzip'
            )
        
        # fill in genotype 
        for hidx, geno_entry in enumerate(genotype_entries):
            if geno_entry not in dosage_row.keys():
                raise ValueError(f'{geno_entry} in genotype_entries are not in dosage_row.')
            self.ARRAY_genotype[hidx, self.n_filled_variants, :] = dosage_row[geno_entry].astype(self.dtype)
            self.n_filled_variants += 1
        
        # fill in variant meta information
        for var_meta_entry, var_meta_list in self.dict_variant_meta.values():
            var_meta_list.append(dosage_row[var_meta_entry])
    
    def _check_equal(self, n1, n2):
        if n1 != n2:
            raise ValueError('Two number are not equal.')
    
    def save(self):
        # run checkers on length before saving
        _check_equal(
            self.nsample, 
            len(self.list_sample)
        )
        for vm in self.dict_variant_meta.keys():
            _check_equal(
                self.nvariant, 
                len(self.dict_variant_meta[vm][1])
            )
        # samples
        self.H5_file.create_dataset(
            'samples', 
            data=self.list_sample, 
            dtype='S25'
        )
        # variants meta
        for kk in self.dict_variant_meta.keys():
            var_meta_entry, var_meta_list = self.dict_variant_meta[kk]
            self.H5_file.create_dataset(
                var_meta_entry, 
                data=var_meta_list, 
                dtype=h5py.string_dtype(encoding='ascii')
            )
        # genotype
        self.H5_file.close()
        
    
        
            
        
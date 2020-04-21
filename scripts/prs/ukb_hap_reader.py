import gc
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()


class UKBhapReader:
    '''
    Reader of BGEN for one chromosome of ukb_hap.
    '''
    def __init__(self, bgen_path, bgen_bgi_path, sample_path, chromosome=None, lazy_load=False):
        if lazy_load is False:
            self.rbgen = importr("rbgen")
            self.bgen_path = bgen_path
            self.bgi_path = bgen_bgi_path
            self.sample_path = sample_path
            self.chromosome = chromosome  # to fix the missing chromosome in ukb_hap_v2
        
    def extract_variant_by_position(self, chrom, start, end, max_entries_per_sample=4):
        '''
        max_entries_per_sample: number of entries per sample in BGEN
        '''
        range_pd = pd.DataFrame(
            {
                'chromosome': [ chrom ],
                'start': [start],
                'end': [end]
            }
        )
#         breakpoint()
        cached_data = self.rbgen.bgen_load(
            self.bgen_path,
            index_filename=self.bgi_path,
            ranges=range_pd, 
            max_entries_per_sample=4
        )
        return cached_data
    
    def _get_varid(self, chrm, pos, a1, a2):
        if self.chromosome is not None and chrm == '':
            chrm = self.chromosome
        return f'{chrm}:{pos}:{a1}:{a2}'
    
    @staticmethod
    def _check_ncol(mat, ncol):
        '''
        check if mat has expected number of columns
        '''
        if mat.shape[1] != ncol:
            raise ValueError(f'ukb_hap does not have {ncol} columns.')
    
    def _hap_to_count(self, hap):
        '''
        expect hap has two columns.
        '''
        self._check_ncol(hap, 2)
        colsum = np.sum(hap, axis=1)
        if (colsum != 1).sum() > 0:
            raise ValueError('some rows have colsum != 1 which is not allowed.')
        return hap[:,1]
    
    @staticmethod
    def _next_pos(curr_pos, max_pos, n_jump):
        next_pos = min(curr_pos + n_jump, max_pos)
        return next_pos
    
    def ukb_hap_to_haplo(self, ukb_hap):
        '''
        expect ukb_hap has 4 columns
        ukb hap encoding law:
            0,1,0,1 -> 1|1
            0,1,1,0 -> 1|0
            1,0,0,1 -> 0|1
            1,0,1,0 -> 0|0
        i.e. the 1st, 2nd columns encode haplotype 1 and 
        the 3rd, 4th columns encode haplotype 2.
        And 0,1->1; 1,0->0. 
        Other combinations are not allowed.
        '''
        self._check_ncol(ukb_hap, 4)
        return self._hap_to_count(ukb_hap[:, :2]), self._hap_to_count(ukb_hap[:, 2:])
    
    def retrieve_from_list(self, chrom, pos, non_effect_allele, effect_allele, 
                           n_var_cached=10, max_entries_per_sample=4):
        '''
        Retrieve generator of variants.
        '''
        niter = 0
        snp_list = pd.DataFrame({
            'chr': chrom,
            'pos': pos,
            'non_effect_allele': non_effect_allele,
            'effect_allele': effect_allele
        })
        # sort by position so that we can retrive by left to right.
        snp_list = snp_list.sort_values('pos').reset_index()
        # set desired variant id 
        set_snp_list = set(
            snp_list.apply(
                lambda x: self._get_varid(x.chr, x.pos, x.non_effect_allele, x.effect_allele),
                axis=1
            ).tolist()
        )
        nsnp = snp_list.shape[0]
        curr_pos_in_snp_list = 0
        next_pos_in_snp_list = 0
        while next_pos_in_snp_list < nsnp:
            # modified from 
            # https://github.com/liangyy/predixcan_prediction/blob/a85d52d89de9fe237a1217b5627c7e8d9f700f7e/bgen/bgen_dosage.py#L80
            next_pos_in_snp_list = self._next_pos(curr_pos_in_snp_list, nsnp, n_var_cached)
            if niter > 0:
                cached_data_struct = cached_data.__sexp__
                del cached_data
                del cached_data_struct
                gc.collect()
            print(snp_list.chr[curr_pos_in_snp_list],snp_list.pos[curr_pos_in_snp_list],snp_list.pos[next_pos_in_snp_list - 1])

            cached_data = self.extract_variant_by_position(
                chrom=snp_list.chr[curr_pos_in_snp_list], 
                start=snp_list.pos[curr_pos_in_snp_list], 
                end=snp_list.pos[next_pos_in_snp_list - 1]
            )
            curr_pos_in_snp_list = next_pos_in_snp_list
            # print(snp_list.chr[curr_pos_in_snp_list],snp_list.pos[curr_pos_in_snp_list],snp_list.pos[next_pos_in_snp_list - 1])
            all_variants = pandas2ri.ri2py(cached_data[0])
            if all_variants.shape[0] == 0:
                return
            all_variants['my_var_id'] = all_variants.apply(
                lambda x: self._get_varid(x.chromosome, x.position, x.allele0, x.allele1),
                axis=1
            )
            all_probs = pandas2ri.ri2py(cached_data[4])
            niter += 1
            for row_idx, (rsid, row) in enumerate(all_variants.iterrows()):
                if row.my_var_id in set_snp_list:
                    dosage_row = row.rename({'chromosome': 'chr'})
                    h1, h2 = self.ukb_hap_to_haplo(all_probs[row_idx, :, :])
                    dosage_row['haplo_dosage_1'] = h1
                    dosage_row['haplo_dosage_2'] = h2
                    yield dosage_row
                    

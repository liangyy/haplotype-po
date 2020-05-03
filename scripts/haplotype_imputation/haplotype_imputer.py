import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()


class HaploImputer:
    
    def __init__(self):
        self.info = '''Imputing parent of origin given PRS and observed phenotype.
        '''
        # self.mode = set('basic_em')
    
    def message(self):
        print(self.info)
    
    def _missing_checker(self, df, desired_values):
        o = df.apply(
            lambda x: self.__is_all_non_missing(x, desired_values), 
            axis=1
        )
        return o
    
    @staticmethod
    def __is_all_non_missing(x, val):
        # breakpoint()
        for xx in x:
            if xx not in val:
                return False
        return True
    
    def __call_rlib_em(self, father, mother, h1, h2, rscript_and_func, return_all=False):
        # all_indivs = h1['individual_id'].tolist()
        df_all = pd.DataFrame({
            'individual_id': h1['individual_id'].tolist()
        })
        
        # filter out individuals with missing observations
        non_miss_in_father = self._missing_checker(
            father.drop('individual_id', axis=1),
            desired_values=[0, 1]
        )
        non_miss_in_mother = self._missing_checker(
            mother.drop('individual_id', axis=1),
            desired_values=[0, 1]
        )
        ff, mm, hh1, hh2 = self._extract_by_rows_binary(
            [father, mother, h1, h2],
            np.logical_and(non_miss_in_father, non_miss_in_mother)
        )

        # drop eid
        fmat, mmat, h1mat, h2mat = self._drop_individual_id(
              [ff, mm, hh1, hh2]
        )
        
        # solve EM
        robjects.r(f'source(\'../../code/{rscript}.R\')'.format(rscript=rscript_and_func[0]))
        em_solver = robjects.globalenv[rscript_and_func[1]]
        
        # output
        out_df = pd.DataFrame({ 'prob_z': pandas2ri.ri2py(out[0]) })
        # breakpoint()
        out_df['individual_id'] = ff['individual_id']
        if return_all is True:
            out_df = pd.merge(
                df_all, out_df, 
                left_on='individual_id',
                right_on='individual_id',
                how='left'
            ).fillna(0.5)
        
        return out_df
    
    def _basic_em(self, father, mother, h1, h2, return_all=False):
        '''
        Only work with individuals has non-missing
        (either 0 or 1) in all columns
        phenotypes in both father and mother.
        Must be called from self.impute. 
        Otherwise the tables may have expected's.
        '''
        return self.__call_rlib_em(
            father, mother, 
            h1, h2, 
            rscript_and_func=('rlib_em.R', 'em_algorithm'), 
            return_all=return_all
        )
    
    def _based_em_deg(self, father, mother, h1, h2, return_all=False):
        '''
        Same as _basic_em
        '''
        return self.__call_rlib_em(
            father, mother, 
            h1, h2, 
            rscript_and_func=('rlib_em_degenerate.R', 'em_algorithm_deg'), 
            return_all=return_all
        )
        
        
        
        
        
        
        
    
    def impute(self, df_father, df_mother, df_h1, df_h2, mode, kwargs={}):
        '''
        wrapper for approaches.
        Individual ID is under individual_id column.
        Work with individuals occur in both of 
        the two PRS tables.
        Work with phenotypes occur in both sets 
        of tables.
        '''
        
        if not hasattr(self, f'_{mode}'):
            raise ValueError(f'{mode} has not been implemented.')
        else:
            impute_method = getattr(self, f'_{mode}')
        
        phenotypes = self._get_common_phenotypes(
            [df_father, df_mother, df_h1, df_h2]
        )
        df_f, df_m, df_1, df_2 = self._extract_by_cols(
            [df_father, df_mother, df_h1, df_h2],
            phenotypes
        )
        
        individuals_prs = self._get_common_individuals(
            [df_1, df_2]
        )
        df_f, df_m, df_1, df_2 = self._extract_and_arrange_rows(
            [df_f, df_m, df_1, df_2],
            individuals_prs
        )
        
        return impute_method(df_f, df_m, df_1, df_2, **kwargs)
    
    @staticmethod
    def _drop_individual_id(dflist):
        o = []
        for df in dflist:
            o.append(df.drop('individual_id', axis=1))
        return o

    @staticmethod
    def _df_to_mat(dflist):
        o = []
        for df in dflist:
            o.append(df.values)
        return o
    
    @staticmethod
    def _extract_by_rows_binary(dflist, binary):
        o = []
        for df in dflist:
            o.append(df[binary].reset_index(drop=True))
        return o
    
    @staticmethod
    def _extract_by_cols(dflist, cols):
        o = []
        for df in dflist:
            o.append(df[cols])
        return o
    
    @staticmethod
    def _extract_and_arrange_rows(dflist, rows):
        ref_df = pd.DataFrame({'individual_id': rows})
        o = []
        for df in dflist:
            o.append(
                pd.merge(
                    ref_df,
                    df,
                    left_on='individual_id',
                    right_on='individual_id',
                    how='left'
                ).reset_index(drop=True)
            )
        return o
        
    def _get_common_phenotypes(self, dflist):
        return self._get_common(dflist, self.__get_columns)
    
    def _get_common_individuals(self, dflist):
        return self._get_common(dflist, self.__get_indivs)
    
    def _get_common(self, dflist, func_to_get_values):
        if len(dflist) < 2:
            raise ValueError('dflist needs more than 1 element.')
        out = set(func_to_get_values(dflist[0]))
        for i in range(1, len(dflist)):
            tmp = set(func_to_get_values(dflist[i]))
            out = tmp & out
        return list(out)
    
    @staticmethod
    def __get_columns(df):
        return df.columns.tolist()
    
    @staticmethod
    def __get_indivs(df):
        return df['individual_id'].tolist()
        
        
        

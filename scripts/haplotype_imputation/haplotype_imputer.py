import pandas as pd
import numpy as np
import torch
from copy import deepcopy
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
        robjects.r('source(\'../../code/{rscript}\')'.format(rscript=rscript_and_func[0]))
        em_solver = robjects.globalenv[rscript_and_func[1]]
        out = em_solver(
            pandas2ri.py2ri(fmat),
            pandas2ri.py2ri(mmat),
            pandas2ri.py2ri(h1mat),
            pandas2ri.py2ri(h2mat)
        ) 

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
    
    @staticmethod
    def _to_torch_tensor(x, device):
        return torch.Tensor(x).to(device)
    
    @staticmethod
    def _check_dim(mat, m, n):
        if mat.shape[0] != m or mat.shape[1] != n:
            raise ValueError(f'Shape does not match the expectation ({m}, {n}).')
    
    def _neg_residual_t_residual_div_2_sigma2(self, r, s2):
        '''
        calculate: - r^t r / (2 * sigma2)
        sigma2: 1 x p
        r: n x p
        '''
        return - torch.sum( self._mat_vec_div_by_row(torch.pow(r, 2) / 2, s2), axis=1 )
    
    @staticmethod
    def _mat_vec_div_by_row(mat, vec):
        '''
        Divide vec to each row of mat: mat_i / vec for all row i
        vec: 1 x p
        '''
        return mat / vec.expand_as(mat)
        
    def _calc_l_cond_z(self, yf, ym, h1, h2, beta, sigma2, z):
        s2f_ = sigma2[0]
        s2m_ = sigma2[1]
        if z == 1:
            yf_ = torch.matmul(h1, beta[0])
            ym_ = torch.matmul(h2, beta[1])
        elif z == 0:
            yf_ = torch.matmul(h2, beta[0])
            ym_ = torch.matmul(h1, beta[1])
        else:
            raise ValueError('Input z can only be 0 or 1.')
            
        rf = yf - yf_  # n x p
        rm = ym - ym_
        
        l = self._neg_residual_t_residual_div_2_sigma2(rf, s2f_) + self._neg_residual_t_residual_div_2_sigma2(rm, s2m_)
        return l
        
    def _calc_l(self, yf, ym, h1, h2, beta, sigma2):
        l0 = self._calc_l_cond_z(yf, ym, h1, h2, beta, sigma2, z = 0)
        l1 = self._calc_l_cond_z(yf, ym, h1, h2, beta, sigma2, z = 1)
        return l0, l1
    
    @staticmethod
    def _calc_gamma(l0, l1):
        max_l = torch.max(l0, l1)
        r1 = l1 - max_l
        r0 = l0 - max_l
        o = torch.exp(r1) / ( torch.exp(r0) + torch.exp(r1) )
        return o
    
    @staticmethod
    def _logsum(a, b):
        '''
        return log(exp(a) + exp(b))
        '''
        m = torch.max(a, b)
        ma = a - m
        mb = b - m
        o = torch.log(torch.exp(ma) + torch.exp(mb)) + m
        return o
    
    @staticmethod
    def _mat_vec_mul_by_col(mat, vec):
        '''
        mat: n x k
        vec: n x 1
        return mat_i / vec for each column i
        '''
        return mat * vec.view(vec.shape[0], -1).expand_as(mat)
    
    @staticmethod
    def _update_beta(N, M, y, pos, beta):
        p = pos.shape[1]
        for pi in range(p):
            mask_ = pos[:, pi]
            y_ = y[:, pi]
            M_ = M[:, mask_]
            N_ = N[mask_, :][:, mask_]
            beta_, _ = torch.solve(
                torch.matmul(M_.T, y_).view(M_.shape[1], -1), 
                N_
            )
            beta[mask_, pi] = beta_[:, 0]
    
    @staticmethod
    def _update_sigma2(n, y, N, beta, sigma2):
        diag_yty = torch.sum(torch.pow(y, 2), axis=0)
        diag_betat_N_beta = torch.einsum('pk,kj,jp->p', beta.T, N, beta)
        sigma2[:] = 1 / n * (diag_yty - diag_betat_N_beta)
    
    def _calc_diff_in_list(self, list_v1, list_v2):
        out = 0
        if len(list_v1) != len(list_v2):
            raise ValueError('Require the input lists to have the same length')
        for i in range(len(list_v1)):
            out += self._calc_diff(list_v1[i], list_v2[i])
        return out
    
    @staticmethod
    def _calc_diff(u, v):
        return torch.sum(torch.pow(u - v, 2))
    
    def _eval_lld(self, l0, l1, sigma2, n):
        '''
        l0: n x 1
        l1: n x 1
        sigma2: 1 x p
        n: scalar
        return log(exp(l0) + exp(l1)) - n / 2 * sum(log(sigma2))
        '''
        l1_plus_l0 = self._logsum(l0, l1)
        o = torch.sum(l1_plus_l0) - n / 2 * torch.sum(torch.log(sigma2[0]) + torch.log(sigma2[1]))
        return o
        
    def _em_otf(self, yf, ym, h1, h2, pos, device='cpu', tol=1e-5, maxiter=100):
        # breakpoint() 
        # add intercept
        h1 = np.concatenate(
            (np.ones((h1.shape[0], 1)), h1),
            axis=1
        )
        h2 = np.concatenate(
            (np.ones((h2.shape[0], 1)), h2),
            axis=1
        )
        
        # push everything to torch Tensor
        h1 = self._to_torch_tensor(h1, device)
        h2 = self._to_torch_tensor(h2, device)
        yf = self._to_torch_tensor(yf, device)
        ym = self._to_torch_tensor(ym, device)
        pos = self._to_torch_tensor(pos, device)

        # add extra rows for intercept in pos 
        pos = torch.cat((torch.ones((1, pos.shape[1])) == 1, pos == 1), axis = 0)
        
        # repeatedly used quantities
        HtH = torch.matmul(h1.T, h1) + torch.matmul(h2.T, h2)
        H = h1 + h2
        
        # dimensions
        p = yf.shape[1]
        n = yf.shape[0]
        k = h1.shape[1]
        self._check_dim(h1, n, k)
        self._check_dim(h2, n, k)
        self._check_dim(yf, n, p)
        self._check_dim(ym, n, p)
        self._check_dim(pos, k, p)
        
        # initilization
        beta = [
            torch.zeros((k, p)),  # father
            torch.zeros((k, p)),  # mother
        ]
        sigma2 = [
            torch.ones((1, p)),  # father
            torch.ones((1, p)),  # mother
        ]
        
        diff = tol + 1
        niter = 0
        lld = []
        # breakpoint()
        while diff > tol and niter < maxiter:
            
            # E step
            l0, l1 = self._calc_l(yf, ym, h1, h2, beta, sigma2)
            # breakpoint()
            lld_curr = self._eval_lld(l0, l1, sigma2, n)
            gamma = self._calc_gamma(l0, l1)
            lld.append(lld_curr)
    
            
            # M step
            
            ## prepare M, N
            d_gamma_h1 = self._mat_vec_mul_by_col(h1, gamma)
            d_ngamma_h2 = self._mat_vec_mul_by_col(h2, (1 - gamma))
            M = d_gamma_h1 + d_ngamma_h2
            N = torch.matmul(h1.T, d_gamma_h1) + torch.matmul(h2.T, d_ngamma_h2)
            tilde_M = H - M
            tilde_N = HtH - N
            
            ## update
            beta_old = deepcopy(beta)
            sigma2_old = deepcopy(sigma2)
            self._update_beta(N, M, yf, pos, beta[0])
            self._update_beta(tilde_N, tilde_M, ym, pos, beta[1])
            self._update_sigma2(n, yf, N, beta[0], sigma2[0])
            self._update_sigma2(n, ym, tilde_N, beta[1], sigma2[1])
            
            # diff
            diff_b = self._calc_diff_in_list(beta_old, beta)
            diff_s = self._calc_diff_in_list(sigma2_old, sigma2)
            
            diff = diff_b + diff_s
            
            # niter
            niter += 1
        
        # last update
        l0, l1 = self._calc_l(yf, ym, h1, h2, beta, sigma2)
        lld_curr = self._eval_lld(l0, l1, sigma2, n)
        gamma = self._calc_gamma(l0, l1)
        lld.append(lld_curr)
        
        return beta, sigma2, gamma, lld
        
    
    def _otf_basic_em(self, father, mother, h1, h2, df_indiv, df_pos, return_all=False):
        '''
        Only work with individuals has non-missing
        (either 0 or 1) in all columns
        phenotypes in both father and mother.
        Must be called from self.impute_otf. 
        Otherwise the tables may not have the expected properties.
        '''
        df_all = pd.DataFrame({
            'individual_id': df_indiv['individual_id'].tolist()
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
        to_keep_ind = np.logical_and(non_miss_in_father, non_miss_in_mother)
        ff, mm = self._extract_by_rows_binary(
            [father, mother],
            to_keep_ind 
        )
        hh1 = h1.T[to_keep_ind, :]
        hh2 = h2.T[to_keep_ind, :]

        # drop eid
        fmat, mmat, posmat = self._drop_individual_id(
              [ff, mm, df_pos]
        )
        
        # solve EM
        # breakpoint()
        # np.save('fmat.npy', fmat.values)
        # np.save('mmat.npy', mmat.values)
        # np.save('hh1.npy', hh1)
        # np.save('hh2.npy', hh2)
        # np.save('posmat.npy', posmat.values)
        beta, sigma2, out, lld = self._em_otf(fmat.values, mmat.values, hh1, hh2, posmat.values)
        
        # output
        out_df = pd.DataFrame({ 'prob_z': out })
        # breakpoint()
        out_df['individual_id'] = ff['individual_id']
        if return_all is True:
            out_df = pd.merge(
                df_all, out_df, 
                left_on='individual_id',
                right_on='individual_id',
                how='left'
            ).fillna(0.5)
        
        return beta, sigma2, out_df, lld
    
    def _basic_em(self, father, mother, h1, h2, return_all=False):
        '''
        Only work with individuals has non-missing
        (either 0 or 1) in all columns
        phenotypes in both father and mother.
        Must be called from self.impute. 
        Otherwise the tables may not have the expected properties.
        '''
        return self.__call_rlib_em(
            father, mother, 
            h1, h2, 
            rscript_and_func=('rlib_em.R', 'em_algorithm'), 
            return_all=return_all
        )
    
    def _basic_em_deg(self, father, mother, h1, h2, return_all=False):
        '''
        Same as _basic_em
        '''
        return self.__call_rlib_em(
            father, mother, 
            h1, h2, 
            rscript_and_func=('rlib_em_degenerate.R', 'em_algorithm_deg'), 
            return_all=return_all
        )
        
        
        
        
        
    def impute_otf(self, df_father, df_mother, h1, h2, indiv_df, pos_df, mode, kwargs={}):
        '''
        wrapper for approaches.
        df_father, df_mother: observed phenotypes.
        h1, h2: haplotypes 
        indiv_list, pos_list: corresponds to the haplotypes.
        snp_dict: for each trait, it gives the SNP list to include.
        '''
        if not hasattr(self, f'_otf_{mode}'):
            raise ValueError(f'{mode} has not been implemented.')
        else:
            impute_method = getattr(self, f'_otf_{mode}')
        
        pos_df['individual_id'] = ''
        phenotypes = self._get_common_phenotypes(
            [df_father, df_mother, pos_df]
        )
        df_f, df_m, df_pos = self._extract_by_cols(
            [df_father, df_mother, pos_df],
            phenotypes
        )
        individuals_hap = indiv_df.individual_id.tolist()
        df_f, df_m, df_indiv = self._extract_and_arrange_rows(
            [df_f, df_m, indiv_df],
            individuals_hap
        )
        return impute_method(df_f, df_m, h1, h2, df_indiv, df_pos, **kwargs)
    
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
        
        
        

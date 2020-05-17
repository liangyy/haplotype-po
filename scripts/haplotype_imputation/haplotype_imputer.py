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
        
    def _em_otf(self, yf, ym, h1, h2, pos, covar=None, device='cpu', tol=1e-5, maxiter=100):
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
        
        # add covariate if there is any
        if covar is not None:
            h1 = np.concatenate(
                (h1, covar),
                axis=1
            )
            h2 = np.concatenate(
                (h2, covar),
                axis=1
            )
            # add column in pos
            # breakpoint()
            pos = pos = np.concatenate( 
                (pos, np.ones((covar.shape[1], pos.shape[1]))), 
                axis=0 
            )
        
        
        # push everything to torch Tensor
        h1 = self._to_torch_tensor(h1, device)
        h2 = self._to_torch_tensor(h2, device)
        yf = self._to_torch_tensor(yf, device)
        ym = self._to_torch_tensor(ym, device)
        pos = self._to_torch_tensor(pos, device)
        
        # breakpoint()

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
            # breakpoint()
        
        # last update
        l0, l1 = self._calc_l(yf, ym, h1, h2, beta, sigma2)
        lld_curr = self._eval_lld(l0, l1, sigma2, n)
        gamma = self._calc_gamma(l0, l1)
        lld.append(lld_curr)
        
        return beta, sigma2, gamma, lld
    
    def _update_beta_per_snp(self, CtC, CtX, XtX, CtY, XtY, pos, beta, beta_c):
        '''
        X: n x k
        Y: n x p
        C: n x Nc
        CtC: Nc x Nc
        CtX: Nc x k
        XtX: k x 1
        CtY: Nc x p
        XtY: k x p
        pos: k x p
        '''
        # check dim
        k = XtX.shape[0]
        Nc = CtC.shape[0]
        p = pos.shape[1]
        self._check_dim(CtC, Nc, Nc)
        self._check_dim(CtX, Nc, k)
        self._check_dim(XtX, k, 1)
        self._check_dim(CtY, Nc, p)
        self._check_dim(XtY, k, p)
        self._check_dim(pos, k, p)
        
        # main
        for pi in range(p):
            mask_ = pos[:, pi]
            CtYp = CtY[:, pi].unsqueeze(axis=1)
            XtXp = XtX[mask_, :]
            CtXp = CtX[:, mask_]
            XtYp = XtY[mask_, :][:, pi].unsqueeze(axis=1)
            beta_, beta_c_ = self.__update_beta_per_snp_one_y(CtC, CtXp, XtXp, CtYp, XtYp)
            # breakpoint()
            beta[mask_, pi] = beta_[:, 0]
            beta_c[:, pi] = beta_c_[:, 0]
    
    @staticmethod
    def __update_beta_per_snp_one_y(CtC, CtX, XtX, CtY, XtY):
        '''
        Equation:
        A = (CtC)^-1 CtY  # Nc x p
        B = (CtC)^-1 CtX  # Nc x k
        D = CtY  # Nc x p
        E = XtY  # k x p
        S = XtX - XtC B  # k x 1 (take diag)
        beta = - S^-1 (Bt D - E)  # k x p
        beta_c = A + B S^-1 (Bt D - E)  # Nc x p
        '''
        A, _ = torch.solve(CtY, CtC)
        B, _ = torch.solve(CtX, CtC)
        D = CtY
        E = XtY
        S = XtX - torch.einsum('nk,nk->k', CtX, B)[:, None]
        # breakpoint()
        # S = torch.unsqueeze(S, axis=1)
        BtD_minus_E = torch.matmul(B.T, D) - E
        beta = - torch.einsum('kl,kp->kp', 1 / S, BtD_minus_E)
        beta_c = A + torch.matmul(B, -beta)
        return beta, beta_c
    
    def _calc_l_per_snp(self, yf, ym, h1, h2, covar_mat, pos, beta, beta_c, sigma2):
        p = pos.shape[1]
        l1 = torch.zeros((yf.shape[0],))
        l0 = torch.zeros((yf.shape[0],))
        for pi in range(p):
            mask_ = pos[:, pi]
            # breakpoint()
            y1_ = yf[:, pi]
            y2_ = ym[:, pi]
            h1_ = h1[:, mask_]
            h2_ = h2[:, mask_]
            beta_f_ = beta[0][mask_, :][:, pi]
            beta_m_ = beta[1][mask_, :][:, pi]
            beta_c_f_ = beta_c[0][:, pi]
            beta_c_m_ = beta_c[1][:, pi]
            sigma2_f_ = sigma2[0][mask_, :][:, pi]
            sigma2_m_ = sigma2[1][mask_, :][:, pi]
            l0_, l1_ = self.__calc_l_per_snp_one_y(y1_, y2_, h1_, h2_, covar_mat, [beta_f_, beta_m_], [beta_c_f_, beta_c_m_], [sigma2_f_, sigma2_m_])
            l1 = l1 + l1_
            l0 = l0 + l0_
        return l0, l1
    
    def __calc_l_per_snp_one_y(self, yf, ym, h1, h2, c, beta, beta_c, sigma2):
        s2f_ = sigma2[0]
        s2m_ = sigma2[1]
        rf1, rf0 = self.__get_residual(yf, yf, h1, h2, c, c, beta[0], beta_c[0])
        rm1, rm0 = self.__get_residual(ym, ym, h2, h1, c, c, beta[1], beta_c[1])
        l1 = self._neg_residual_t_residual_div_2_sigma2_per_snp(rf1, s2f_) + self._neg_residual_t_residual_div_2_sigma2_per_snp(rm1, s2m_)
        l0 = self._neg_residual_t_residual_div_2_sigma2_per_snp(rf0, s2f_) + self._neg_residual_t_residual_div_2_sigma2_per_snp(rm0, s2m_)
        # breakpoint()
        return l0, l1
    
    def _neg_residual_t_residual_div_2_sigma2_per_snp(self, res, s2):
        '''
        res: n x k
        s2: k
        return: - res ** 2 / 2 / s2
        '''
        # breakpoint()
        o = - torch.pow(res, 2) / 2 / s2[None, :]
        return o.sum(axis=1)
    
    def _update_sigma2_per_snp(self, n, y1, y2, h1, h2, c1, c2, pos, beta, beta_c, sigma2):
        p = pos.shape[1]
        for pi in range(p):
            mask_ = pos[:, pi]
            y1_ = y1[:, pi]
            y2_ = y2[:, pi]
            h1_ = h1[:, mask_]
            h2_ = h2[:, mask_]
            beta_ = beta[:, pi][mask_]
            beta_c_ = beta_c[:, pi]
            sigma2_ = self.__update_sigma2_per_snp_one_y(n, y1_, y2_, h1_, h2_, c1, c2, beta_, beta_c_)
            sigma2[mask_, pi] = sigma2_
    
    @staticmethod
    def __get_residual(y1, y2, x1, x2, c1, c2, beta, beta_c):
        '''
        Equation
        
        cb1 = c1 beta_c
        cb2 = c2 beta_c
        
        for each column of h1, and h2
        
        r1 = y1 - h1 beta - cb1
        r2 = y2 - h2 beta - cb2
        '''
        yp1 = torch.einsum('nk,k->nk', x1, beta)
        yp2 = torch.einsum('nk,k->nk', x2, beta)
        cb1 = torch.matmul(c1, beta_c)
        cb2 = torch.matmul(c2, beta_c)
        
        r1 = y1[:, None] - yp1 - cb1[:, None]
        r2 = y2[:, None] - yp2 - cb2[:, None]
        return r1, r2
    
    def __update_sigma2_per_snp_one_y(self, n, y1, y2, h1, h2, c1, c2, beta, beta_c):
        '''
        get residual from __get_residual
        rsq = r1.T r1 + r2.T r2
        sigma2 = rsq / n
        '''
        r1, r2 = self.__get_residual(y1, y2, h1, h2, c1, c2, beta, beta_c)
        rsq = torch.pow(r1, 2) + torch.pow(r2, 2)
        sigma2 = rsq.sum(axis=0) / n
        
        return sigma2
    
    def _eval_lld_per_snp(self, l0, l1, sigma2, n):
        return self._eval_lld(l0, l1, sigma2, n)
    
    def _calc_gamma_per_snp(self, l0, l1):
        return self._calc_gamma(l0, l1)
        
    def _update_beta_and_sigma2_per_snp(self, beta, beta_c, sigma2, h1, h2, yf, ym, gamma, pos, covar_mat, HtH, CtH, CtC, CtYf, CtYm, n):
        ## prepare CtC, CtX, XtX, CtY, XtY
        ## and update beta 
        d_sqrt_gamma_h1 = self._mat_vec_mul_by_col(h1, torch.sqrt(gamma))
        d_sqrt_ngamma_h2 = self._mat_vec_mul_by_col(h2, torch.sqrt(1 - gamma))
        d_sqrt_ngamma_h1 = self._mat_vec_mul_by_col(h1, torch.sqrt(1 - gamma))
        d_sqrt_gamma_h2 = self._mat_vec_mul_by_col(h2, torch.sqrt(gamma))
        d_sqrt_gamma_covar = self._mat_vec_mul_by_col(covar_mat, torch.sqrt(gamma))
        d_sqrt_ngamma_covar = self._mat_vec_mul_by_col(covar_mat, torch.sqrt(1 - gamma))
        ### for father
        d_sqrt_gamma_yf = self._mat_vec_mul_by_col(yf, torch.sqrt(gamma))
        d_sqrt_ngamma_yf = self._mat_vec_mul_by_col(yf, torch.sqrt(1 - gamma))
        CtX = torch.matmul(d_sqrt_gamma_covar.T, d_sqrt_gamma_h1) + torch.matmul(d_sqrt_ngamma_covar.T, d_sqrt_ngamma_h2)
        XtX = torch.einsum('ij,ij->j', d_sqrt_gamma_h1, d_sqrt_gamma_h1) + torch.einsum('ij,ij->j', d_sqrt_ngamma_h2, d_sqrt_ngamma_h2)
        XtX = torch.unsqueeze(XtX, axis=1)
        XtYf = torch.matmul(d_sqrt_gamma_h1.T, d_sqrt_gamma_yf) + torch.matmul(d_sqrt_ngamma_h2.T, d_sqrt_ngamma_yf)
        self._update_beta_per_snp(CtC, CtX, XtX, CtYf, XtYf, pos, beta[0], beta_c[0])
        ### for mother
        d_sqrt_gamma_ym = self._mat_vec_mul_by_col(ym, torch.sqrt(gamma))
        d_sqrt_ngamma_ym = self._mat_vec_mul_by_col(ym, torch.sqrt(1 - gamma))
        CtX = CtH - CtX
        XtX = HtH - XtX
        XtYm = torch.matmul(d_sqrt_ngamma_h1.T, d_sqrt_ngamma_ym) + torch.matmul(d_sqrt_gamma_h2.T, d_sqrt_gamma_ym)
        # XtYm = HtYm - XtYm
        self._update_beta_per_snp(CtC, CtX, XtX, CtYm, XtYm, pos, beta[1], beta_c[1])
        
        ## update sigma
        self._update_sigma2_per_snp(n, d_sqrt_gamma_yf, d_sqrt_ngamma_yf, d_sqrt_gamma_h1, d_sqrt_ngamma_h2, d_sqrt_gamma_covar, d_sqrt_ngamma_covar, pos, beta[0], beta_c[0], sigma2[0])
        self._update_sigma2_per_snp(n, d_sqrt_gamma_ym, d_sqrt_ngamma_ym, d_sqrt_gamma_h2, d_sqrt_ngamma_h1, d_sqrt_gamma_covar, d_sqrt_ngamma_covar, pos, beta[1], beta_c[1], sigma2[1])
    
    
    def _em_otf_per_snp(self, yf, ym, h1, h2, pos, covar=None, device='cpu', tol=1e-5, maxiter=100):
        # breakpoint() 
        # add intercept as part of covariate matrix
        covar_mat = np.ones((h1.shape[0], 1))
        
        # add covariate if there is any
        if covar is not None:
            covar_mat = np.concatenate(
                (covar_mat, covar),
                axis=1
            )
        
        
        # push everything to torch Tensor
        h1 = self._to_torch_tensor(h1, device)
        h2 = self._to_torch_tensor(h2, device)
        yf = self._to_torch_tensor(yf, device)
        ym = self._to_torch_tensor(ym, device)
        pos = self._to_torch_tensor(pos, device) == 1
        covar_mat = self._to_torch_tensor(covar_mat, device)
        
        
        # repeatedly used terms
        HtH = torch.einsum('nk,nk->k', h1, h1) + torch.einsum('nk,nk->k', h2, h2)
        HtH = torch.unsqueeze(HtH, axis=1)
        # HtYm = torch.matmul(h1.T, ym) + torch.matmul(h2.T, ym)
        CtH = torch.matmul(covar_mat.T, h1) + torch.matmul(covar_mat.T, h2)
        CtC = torch.matmul(covar_mat.T, covar_mat)
        CtYf = torch.matmul(covar_mat.T, yf)
        CtYm = torch.matmul(covar_mat.T, ym)
        
        
        # dimensions
        p = yf.shape[1]
        n = yf.shape[0]
        k = h1.shape[1]
        ncovar = covar_mat.shape[1]
        self._check_dim(h1, n, k)
        self._check_dim(h2, n, k)
        self._check_dim(yf, n, p)
        self._check_dim(ym, n, p)
        self._check_dim(pos, k, p)
        self._check_dim(covar_mat, n, ncovar)
        
        # initilization
        beta = [
            torch.zeros((k, p)),  # father
            torch.zeros((k, p)),  # mother
        ]
        beta_c = [
            torch.zeros((ncovar, p)),  # father
            torch.zeros((ncovar, p)),  # mother
        ]
        sigma2 = [
            torch.ones((k, p)),  # father
            torch.ones((k, p)),  # mother
        ]
        
        diff = tol + 1
        niter = 0
        lld = []
        # breakpoint()
        while diff > tol and niter < maxiter:
            
            # E step
            l0, l1 = self._calc_l_per_snp(yf, ym, h1, h2, covar_mat, pos, beta, beta_c, sigma2)
            # breakpoint()
            lld_curr = self._eval_lld_per_snp(l0, l1, sigma2, n)
            gamma = self._calc_gamma_per_snp(l0, l1)
            lld.append(lld_curr)
            # breakpoint()
    
            
            # M step
            
            ## cache variables
            beta_old = deepcopy(beta)
            sigma2_old = deepcopy(sigma2)
            
            ## update beta and sigma2
            # breakpoint()
            self._update_beta_and_sigma2_per_snp(
                beta, beta_c, sigma2, 
                h1, h2, yf, ym, gamma, pos, covar_mat, 
                HtH, CtH, CtC, CtYf, CtYm, n
            )
            
            # diff
            diff_b = self._calc_diff_in_list(beta_old, beta)
            diff_s = self._calc_diff_in_list(sigma2_old, sigma2)
            
            diff = diff_b + diff_s
            
            # niter
            niter += 1
            # breakpoint()
        
        # last update
        l0, l1 = self._calc_l_per_snp(yf, ym, h1, h2, covar_mat, pos, beta, beta_c, sigma2)
        lld_curr = self._eval_lld_per_snp(l0, l1, sigma2, n)
        gamma = self._calc_gamma_per_snp(l0, l1)
        lld.append(lld_curr)
        
        return beta, beta_c, sigma2, gamma, lld
        
    
    def _otf_basic_em(self, father, mother, h1, h2, df_indiv, df_pos, df_covar=None, return_all=False):
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
        
        if df_covar is None:
            ff, mm = self._extract_by_rows_binary(
                [father, mother],
                to_keep_ind 
            )
            cc = None
        else:
            ff, mm, cc = self._extract_by_rows_binary(
                [father, mother, df_covar],
                to_keep_ind 
            )
        hh1 = h1.T[to_keep_ind, :]
        hh2 = h2.T[to_keep_ind, :]

        # drop eid
        if cc is None:
            fmat, mmat, posmat = self._drop_individual_id(
                  [ff, mm, df_pos]
            )
            cmat = None
        else:
            fmat, mmat, cmat, posmat = self._drop_individual_id(
                  [ff, mm, cc, df_pos]
            )
        
        # solve EM
        # breakpoint()
        # np.save('fmat.npy', fmat.values)
        # np.save('mmat.npy', mmat.values)
        # np.save('hh1.npy', hh1)
        # np.save('hh2.npy', hh2)
        # np.save('posmat.npy', posmat.values)
        # beta, sigma2, out, lld = self._em_otf(fmat.values, mmat.values, hh1, hh2, posmat.values)
        beta, sigma2, out, lld = self._em_otf(fmat.values, mmat.values, hh1, hh2, posmat.values, covar=cmat)
        
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
        
        
        
        
        
    def impute_otf(self, df_father, df_mother, h1, h2, indiv_df, pos_df, mode, df_covar=None, kwargs={}):
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
        
        if df_covar is None:
            df_f, df_m, df_indiv = self._extract_and_arrange_rows(
                [df_f, df_m, indiv_df],
                individuals_hap
            )
            df_c = None
        else:
            df_f, df_m, df_c, df_indiv = self._extract_and_arrange_rows(
                [df_f, df_m, df_covar, indiv_df],
                individuals_hap
            )
            
        return impute_method(df_f, df_m, h1, h2, df_indiv, df_pos, df_covar=df_c, **kwargs)
    
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
        
        
        

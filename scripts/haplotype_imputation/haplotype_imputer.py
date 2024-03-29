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
    
    def _nan_checker(self, df):
        o = df.apply(
            lambda x: self.__is_all_non_nan(x),
            axis=1
        )
        return o

    @staticmethod
    def __is_all_non_nan(x):
        return x.isnull().sum() == 0

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
    def _mat_vec_add_by_row(mat, vec):
        '''
        Add vec to each row of mat: mat_i + vec for all row i
        vec: 1 x p
        '''
        return mat + vec.expand_as(mat)
    
    @staticmethod
    def _mat_vec_div_by_row(mat, vec):
        '''
        Divide vec to each row of mat: mat_i / vec for all row i
        vec: 1 x p
        '''
        return mat / vec.expand_as(mat)
    
    @staticmethod
    def _mat_vec_mul_by_row(mat, vec):
        '''
        Add vec to each row of mat: mat_i + vec for all row i
        vec: 1 x p
        '''
        return mat * vec.expand_as(mat)
    
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
    
    def _eval_lld(self, l0, l1, sigma2, n, prior_z=0.5):
        '''
        l0: n x 1
        l1: n x 1
        sigma2: 1 x p
        n: scalar
        return log(exp(l0) + exp(l1)) - n / 2 * sum(log(sigma2))
        '''
        # FIXME: not sure what's the best way to have prior_z integrated ..
        # prior_z_t = torch.ones((1)) * prior_z
        # l1_plus_l0 = self._logsum(l0 + torch.log(1 - prior_z_t), l1 + torch.log(prior_z_t))
        # END
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
            beta_c[:, mask_, pi] = beta_c_[:, :, 0]
    
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
        beta_c = A + B S^-1 (Bt D - E)  # Nc x k x p
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
        beta_c = A[:, None, :] + torch.einsum('nk,kp->nkp', B, -beta)
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
            beta_c_f_ = beta_c[0][:, mask_, :][:, :, pi]
            beta_c_m_ = beta_c[1][:, mask_, :][:, :, pi]
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
            beta_c_ = beta_c[:, :, pi][:, mask_]
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
        cb1 = torch.einsum('nc,ck->nk', c1, beta_c)
        cb2 = torch.einsum('nc,ck->nk', c2, beta_c)
        
        r1 = y1[:, None] - yp1 - cb1
        r2 = y2[:, None] - yp2 - cb2
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
            torch.zeros((ncovar, k, p)),  # father
            torch.zeros((ncovar, k, p)),  # mother
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
       
        # beta[0] = torch.cat((beta_c[0], beta[0]), axis=0)
        # beta[1] = torch.cat((beta_c[1], beta[1]), axis=0)

        return beta, beta_c, sigma2, gamma, lld
    
    def __call_otf_em(self, father, mother, h1, h2, df_indiv, df_pos, em_func, df_covar=None, return_all=False, debug_cache=None): 
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
        
        # remove SNPs with constant dosage in any haplotype
        non_const_dos_ind = self._get_constant_snp([hh1, hh2])
        hh1 = hh1[:, non_const_dos_ind]
        hh2 = hh2[:, non_const_dos_ind]
        posmat = posmat[non_const_dos_ind].reset_index(drop=True)
        
        # solve EM
        if debug_cache is not None:
        # breakpoint()
            np.save(f'{debug_cache}_fmat.npy', fmat.values)
            np.save(f'{debug_cache}_mmat.npy', mmat.values)
            np.save(f'{debug_cache}_hh1.npy', hh1)
            np.save(f'{debug_cache}_hh2.npy', hh2)
            np.save(f'{debug_cache}_posmat.npy', posmat.values)
            np.save(f'{debug_cache}_cmat.npy', cmat)
        # beta, sigma2, out, lld = self._em_otf(fmat.values, mmat.values, hh1, hh2, posmat.values)
        return em_func(fmat.values, mmat.values, hh1, hh2, posmat.values, covar=cmat), ff['individual_id']
    
    def _otf_per_snp_em(self, father, mother, h1, h2, df_indiv, df_pos, df_covar=None, return_all=False, debug_cache=None):
        '''
        Only work with individuals has non-missing
        (either 0 or 1) in all columns
        phenotypes in both father and mother.
        Must be called from self.impute_otf. 
        Otherwise the tables may not have the expected properties.
        '''
        (beta, beta_c, sigma2, out, lld), indiv_id = self.__call_otf_em(father, mother, h1, h2, df_indiv, df_pos, em_func=self._em_otf_per_snp, df_covar=df_covar, return_all=return_all, debug_cache=debug_cache)
        # beta[0] = torch.cat((beta_c[0], beta[0]), axis=0)
        # beta[1] = torch.cat((beta_c[1], beta[1]), axis=0)
        
        # output
        out_df = pd.DataFrame({ 'prob_z': out })
        # breakpoint()
        out_df['individual_id'] = indiv_id
        if return_all is True:
            df_all = pd.DataFrame({
                'individual_id': df_indiv['individual_id'].tolist()
            })
            out_df = pd.merge(
                df_all, out_df, 
                left_on='individual_id',
                right_on='individual_id',
                how='left'
            ).fillna(0.5)
        
        return (beta, beta_c), sigma2, out_df, lld
        
    def _otf_basic_em(self, father, mother, h1, h2, df_indiv, df_pos, df_covar=None, return_all=False, debug_cache=None):
        '''
        Only work with individuals has non-missing
        (either 0 or 1) in all columns
        phenotypes in both father and mother.
        Must be called from self.impute_otf. 
        Otherwise the tables may not have the expected properties.
        '''
        (beta, sigma2, out, lld), indiv_id = self.__call_otf_em(father, mother, h1, h2, df_indiv, df_pos, em_func=self._em_otf, df_covar=df_covar, return_all=return_all, debug_cache=debug_cache)
        
        # output
        out_df = pd.DataFrame({ 'prob_z': out })
        # breakpoint()
        out_df['individual_id'] = indiv_id
        if return_all is True:
            df_all = pd.DataFrame({
                'individual_id': df_indiv['individual_id'].tolist()
            })
            out_df = pd.merge(
                df_all, out_df, 
                left_on='individual_id',
                right_on='individual_id',
                how='left'
            ).fillna(0.5)
        
        return beta, sigma2, out_df, lld
    
    def __call_em_py(self, father, mother, h1, h2, em_func, df_covar=None, debug_cache=None):
        df_all = pd.DataFrame({
            'individual_id': h1['individual_id'].tolist()
        })
        
        # filter out individuals with missing observations
        non_miss_in_father = self._nan_checker(
            father.drop('individual_id', axis=1)
        )
        non_miss_in_mother = self._nan_checker(
            mother.drop('individual_id', axis=1)
        )
        # breakpoint() 
        to_keep_ind = np.logical_and(non_miss_in_father, non_miss_in_mother)
        
        if df_covar is None:
            ff, mm, hh1, hh2 = self._extract_by_rows_binary(
                [father, mother, h1, h2],
                to_keep_ind 
            )
            cc = None
        else:
            ff, mm, hh1, hh2, cc = self._extract_by_rows_binary(
                [father, mother, h1, h2, df_covar],
                to_keep_ind 
            )
        
        # drop eid
        if cc is None:
            fmat, mmat, h1mat, h2mat = self._drop_individual_id(
                  [ff, mm, hh1, hh2]
            )
            cmat = None
        else:
            fmat, mmat, h1mat, h2mat, cmat = self._drop_individual_id(
                  [ff, mm, hh1, hh2, cc]
            )
            
        # solve EM
        if debug_cache is not None:
        # breakpoint()
            np.save(f'{debug_cache}_fmat.npy', fmat.values)
            np.save(f'{debug_cache}_mmat.npy', mmat.values)
            np.save(f'{debug_cache}_h1mat.npy', h1mat.values)
            np.save(f'{debug_cache}_h2mat.npy', h2mat.values)
            np.save(f'{debug_cache}_cmat.npy', cmat.values)
       
        # breakpoint()
        return em_func(fmat.values, mmat.values, h1mat.values, h2mat.values, covar=cmat), ff['individual_id']
    
    def _avar_update_sigma(self, y, g_1, g_2, beta, cc, omega):
        r1 = self._avar_get_residual(y, g_1, cc, beta)
        r2 = self._avar_get_residual(y, g_2, cc, beta)
        r_tilde = self._mat_vec_mul_by_col(torch.cat((r1, r2), axis=0), torch.sqrt(omega))
        denom = torch.sum(omega)
        sigma2 = torch.sum(torch.pow(r_tilde, 2), axis=0) / denom
        return sigma2
    
    def _avar_update_beta(self, y, g_1, g_2, cc, omega, non_negative):
        # update all phenotypes simultaneously
        # Equation:
        #   Y: n x p
        #   X: n x p
        #   C: n x Nc
        #   A = (CtC)^-1 CtY  # Nc x p
        #   B = (CtC)^-1 CtX  # Nc x p
        #   D = CtY  # Nc x p
        #   E = XtY  # p x 1 (take diag)
        #   S = XtX - XtC B  # p x 1 (take diag)  
        #   BtD = Bt D  # p x 1 (take diag)
        #   beta = - S^-1 (BtD - E)  # p x 1
        #   beta_c = A - B beta (by row)  # Nc x p
        CtC = torch.matmul(cc.T, cc)
        CtY = torch.matmul(cc.T, y)
        Y = self._mat_vec_mul_by_col(torch.cat((y, y), axis=0), torch.sqrt(omega))
        X = self._mat_vec_mul_by_col(torch.cat((g_1, g_2), axis=0), torch.sqrt(omega))
        C = self._mat_vec_mul_by_col(torch.cat((cc, cc), axis=0), torch.sqrt(omega))
        CtX = torch.matmul(C.T, X)
        A, _ = torch.solve(CtY, CtC)
        B, _ = torch.solve(CtX, CtC)
        D = CtY
        E = torch.sum(X * Y, axis=0)
        # breakpoint()
        S = torch.sum(X * X, axis=0) - torch.sum(CtX * B, axis=0)
        BtD = torch.sum(B * D, axis=0)
        beta = - (BtD - E) / S
        # the extra constraint on beta so that beta is non negative
        if non_negative is True:
            beta[beta < 0] = 0
        
        beta_c = A - self._mat_vec_mul_by_row(B, beta)
        # breakpoint()
        return torch.cat((torch.unsqueeze(beta, axis=0), beta_c), axis=0)
    
    def _avar_get_lld(self, l1, l0):
        return torch.sum(self._logsum(-l1, -l0))
    
    def _avar_get_residual(self, y, g, cc, beta):
        yg = self._mat_vec_mul_by_row(g, beta[0, :])
        ycovar = torch.matmul(cc, beta[1:, :])
        return y - yg - ycovar
    
    def _avar_nlog_prob_y_given_z(self, yf, ym, gf, gm, covar, beta, sigma2):
        res_f = self._avar_get_residual(yf, gf, covar, beta[0])
        res_m = self._avar_get_residual(ym, gm, covar, beta[1])
        ratio_f = self._mat_vec_div_by_row(res_f ** 2, 2 * sigma2[0])
        ratio_m = self._mat_vec_div_by_row(res_m ** 2, 2 * sigma2[1])
        lf_n_by_p = self._mat_vec_add_by_row(ratio_f, 1 / 2 * torch.log(sigma2[0]))
        lm_n_by_p = self._mat_vec_add_by_row(ratio_m, 1 / 2 * torch.log(sigma2[1]))
        return torch.sum(lf_n_by_p, axis=1) + torch.sum(lm_n_by_p, axis=1)
    
    def _em_py(self, yf, ym, h1, h2, covar=None, device='cpu', tol=1e-5, maxiter=100, non_negative=True):
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
        covar_mat = self._to_torch_tensor(covar_mat, device)
        
        # dimensions
        p = yf.shape[1]
        n = yf.shape[0]
        k = h1.shape[1]
        ncovar = covar_mat.shape[1]
        self._check_dim(h1, n, k)
        self._check_dim(h2, n, k)
        self._check_dim(yf, n, p)
        self._check_dim(ym, n, p)
        self._check_dim(covar_mat, n, ncovar) 
        
        # initilization
        ## beta: PRS coefficient goes first
        beta = [
            torch.zeros((1 + ncovar, p)),  # father
            torch.zeros((1 + ncovar, p)),  # mother
        ]
        sigma2 = [
            torch.ones((p,)),  # father
            torch.ones((p,)),  # mother
        ]
        
        diff = tol + 1
        niter = 0
        lld = []
        
        while diff > tol and niter < maxiter:
            
            # E step
            l1 = self._avar_nlog_prob_y_given_z(
                yf, ym, 
                h1, h2, covar_mat,
                beta, sigma2
            )
            l0 = self._avar_nlog_prob_y_given_z(
                yf, ym, 
                h2, h1, covar_mat,
                beta, sigma2
            )
            # breakpoint()
            lld.append(self._avar_get_lld(l1, l0))
            gamma = 1 / (1 + torch.exp(l1 - l0))

            # M step
            beta_old = deepcopy(beta)
            sigma2_old = deepcopy(sigma2)
            beta[0] = self._avar_update_beta(
                yf, h1, h2, covar_mat, 
                omega=torch.cat((gamma, 1 - gamma)), 
                non_negative=non_negative
            )
            beta[1] = self._avar_update_beta(
                ym, h2, h1, covar_mat, 
                omega = torch.cat((gamma, 1 - gamma)),
                non_negative=non_negative
            )
            sigma2[0] = self._avar_update_sigma(
                yf, h1, h2, beta[0], covar_mat, 
                omega=torch.cat((gamma, 1 - gamma))
            )
            sigma2[1] = self._avar_update_sigma(
                ym, h2, h1, beta[1], covar_mat, 
                omega=torch.cat((gamma, 1 - gamma))
            )
            
            # calc diff
            diff_b = self._calc_diff_in_list(beta_old, beta)
            diff_s = self._calc_diff_in_list(sigma2_old, sigma2)
            
            diff = diff_b + diff_s
            
            # niter
            niter += 1
            
        return beta, sigma2, gamma, lld
        
    def _basic_em_py(self, father, mother, h1, h2, covar=None, debug_cache=None):
        '''
        Only work with individuals has non-missing
        (either 0 or 1) in all columns
        phenotypes in both father and mother.
        Must be called from self.impute. 
        Otherwise the tables may not have the expected properties.
        '''
        (beta, sigma2, out, lld), indiv_id = self.__call_em_py(father, mother, h1, h2, em_func=self._em_py, df_covar=covar, debug_cache=debug_cache)
        # output
        out_df = pd.DataFrame({ 'prob_z': out })
        # breakpoint()
        out_df['individual_id'] = indiv_id
        
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
    
    def _check_if_trait_in_list(self, traits, mylist):
        for trait in traits:
            if trait not in mylist:
                return False
        return True
    
    def impute_preload_genotype(self, h1, h2, indiv_df, pos_df, phenotypes, individual_ids, output_prefix):
        # breakpoint()
        if not self._check_if_trait_in_list(pos_df.columns.tolist(), phenotypes):
            print('Some phenotypes in corresponding pre-loaded NPY are missed in the input variant dict.')
            return None

        posmat = self._extract_by_cols(
            [pos_df],
            phenotypes
        )[0]
        
        rearrange_idx = self._get_match_idx(indiv_df, individual_ids)
        hh1 = h1.T[rearrange_idx, :]
        hh2 = h2.T[rearrange_idx, :]

        # remove SNPs with constant dosage in any haplotype
        non_const_dos_ind = self._get_constant_snp([hh1, hh2])
        hh1 = hh1[:, non_const_dos_ind]
        hh2 = hh2[:, non_const_dos_ind]
        # breakpoint()
        posmat = posmat[non_const_dos_ind].reset_index(drop=True)
        
        # save
        np.save(f'{output_prefix}.hh1.npy', np.int8(hh1))
        np.save(f'{output_prefix}.hh2.npy', np.int8(hh2))
        np.save(f'{output_prefix}.posmat.npy', posmat.values)
        
        
        
    def impute_preload_pheno_and_covar(self, df_father, df_mother, indiv_df, output_prefix, df_covar=None):
        
        # stage 1
        phenotypes = self._get_common_phenotypes(
            [df_father, df_mother]
        )
        df_f, df_m = self._extract_by_cols(
            [df_father, df_mother],
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
        
        # stage 2
        non_miss_in_father = self._missing_checker(
            df_f.drop('individual_id', axis=1),
            desired_values=[0, 1]
        )
        non_miss_in_mother = self._missing_checker(
            df_m.drop('individual_id', axis=1),
            desired_values=[0, 1]
        )
            
        to_keep_ind = np.logical_and(non_miss_in_father, non_miss_in_mother)
        
        if df_c is None:
            ff, mm = self._extract_by_rows_binary(
                [df_f, df_m],
                to_keep_ind 
            )
            cc = None
        else:
            ff, mm, cc = self._extract_by_rows_binary(
                [df_f, df_m, df_c],
                to_keep_ind 
            )

        # drop eid
        if cc is None:
            fmat, mmat = self._drop_individual_id(
                  [ff, mm]
            )
            cmat = None
        else:
            fmat, mmat, cmat = self._drop_individual_id(
                  [ff, mm, cc]
            )
            np.save(f'{output_prefix}.cmat.npy', cmat.values)
        np.save(f'{output_prefix}.fmat.npy', fmat.values)
        np.save(f'{output_prefix}.mmat.npy', mmat.values)
        np.save(f'{output_prefix}.individual_id.npy', ff['individual_id'].values)
        np.save(f'{output_prefix}.phenotype.npy', fmat.columns.tolist())    

            
    def impute_otf_multi_chr(self, yf, ym, pos_dict, chroms, genotype1_pattern, genotype2_pattern, cmat=None, device='cpu', tol=1e-5, maxiter=100):
        # initialize
        ## add intercept and covariates
        covar = np.ones((yf.shape[0], 1))
        pos_covar = np.ones((1, yf.shape[1]))
        if cmat is not None:
            covar = np.concatenate(
                (covar, cmat),
                axis=1
            )
            pos_covar = np.concatenate( 
                (pos_covar, np.ones((cmat.shape[1], yf.shape[1]))), 
                axis=0 
            )
        
        ## push everything to torch Tensor
        yf = self._to_torch_tensor(yf, device)
        ym = self._to_torch_tensor(ym, device)
        for chrom in chroms:
            pos_dict[chrom] = self._to_torch_tensor(pos_dict[chrom], device)
        if cmat is not None:
            covar = self._to_torch_tensor(covar, device)
            pos_covar = self._to_torch_tensor(pos_covar, device)
        
        ## getting and checking dimensions
        p = yf.shape[1]
        n = yf.shape[0]
        k_dict = { chrom: pos_dict[chrom].shape[0] for chrom in chroms }
        kcovar = pos_covar.shape[0]
        self._check_dim(yf, n, p)
        self._check_dim(ym, n, p)
        self._check_dim(pos_covar, kcovar, p)
        for chrom in chroms:
            self._check_dim(pos_dict[chrom], k_dict[chrom], p)
        
        ## combine pos
        for chrom in chroms:
            pos_dict[chrom] = torch.cat((pos_dict[chrom] == 1, pos_covar == 1), axis=0)
        
        ## parameters
        beta = { 
            chrom: [
                torch.zeros((k_dict[chrom] + kcovar, p)),  # father
                torch.zeros((k_dict[chrom] + kcovar, p)),  # mother
            ] for chrom in chroms
        }
        sigma2 = { 
            chrom: [
                torch.ones((1, p)),  # father
                torch.ones((1, p)),  # mother
            ] for chrom in chroms
        }
        # breakpoint() 
        ## initialize leave-one-chromsome-out residual (LOCOR)
        ## importantly, we initialized betas as zeros so that LOCOR equals to observed phenotype
        locorf = yf
        locorm = ym 
        
        # END of initilization
        
        
        diff = tol + 1
        niter = 0
        gamma_list = { chrom: None for chrom in chroms }
        # here we maintain per chromosome lld since the all-chromosome lld is not tractable
        # FIXME: will it always decrease?
        lld = { chrom: [] for chrom in chroms }
        
        while diff > tol and niter < maxiter:
            
            diff_b = 0
            diff_s = 0
            for chrom in chroms:
                # load haplotype
                h1 = self._load_haplotype_from_preloaded(genotype1_pattern.format(chr_num=chrom), device)
                h2 = self._load_haplotype_from_preloaded(genotype2_pattern.format(chr_num=chrom), device)
                h1 = torch.cat((h1, covar), axis=1)
                h2 = torch.cat((h2, covar), axis=1)
                # breakpoint()

                # add back genetic effect of current chromosome
                locorf = locorf + self._get_avg_genetic_effect(h1[:, :k_dict[chrom]], h2[:, :k_dict[chrom]], beta[chrom][0][:k_dict[chrom], :])
                locorm = locorm + self._get_avg_genetic_effect(h1[:, :k_dict[chrom]], h2[:, :k_dict[chrom]], beta[chrom][1][:k_dict[chrom], :])
                
                # E step
                l0, l1 = self._calc_l(locorf, locorm, h1, h2, beta[chrom], sigma2[chrom])
                # breakpoint()
                lld_curr = self._eval_lld(l0, l1, sigma2[chrom], n)
                gamma_list[chrom] = self._calc_gamma(l0, l1)
                gamma = gamma_list[chrom]
                lld[chrom].append(lld_curr)
                
                # M step
                
                ## to save memory we don't keep these variables but do inline calculation
                # HtH = torch.matmul(h1.T, h1) + torch.matmul(h2.T, h2)
                # H = h1 + h2
                
                ## prepare M, N
                d_gamma_h1 = self._mat_vec_mul_by_col(h1, gamma)
                d_ngamma_h2 = self._mat_vec_mul_by_col(h2, (1 - gamma))
                M = d_gamma_h1 + d_ngamma_h2
                N = torch.matmul(h1.T, d_gamma_h1) + torch.matmul(h2.T, d_ngamma_h2)
                tilde_M = h1 + h2 - M
                tilde_N = torch.matmul(h1.T, h1) + torch.matmul(h2.T, h2) - N
                
                ## update
                beta_old = deepcopy(beta[chrom])
                sigma2_old = deepcopy(sigma2[chrom])
                self._update_beta(N, M, locorf, pos_dict[chrom], beta[chrom][0])
                self._update_beta(tilde_N, tilde_M, locorm, pos_dict[chrom], beta[chrom][1])
                self._update_sigma2(n, locorf, N, beta[chrom][0], sigma2[chrom][0])
                self._update_sigma2(n, locorm, tilde_N, beta[chrom][1], sigma2[chrom][1])
                # breakpoint()

                # subtract genetic effect of current chromosome
                locorf = locorf - self._get_avg_genetic_effect(h1[:, :k_dict[chrom]], h2[:, :k_dict[chrom]], beta[chrom][0][:k_dict[chrom], :])
                locorm = locorm - self._get_avg_genetic_effect(h1[:, :k_dict[chrom]], h2[:, :k_dict[chrom]], beta[chrom][1][:k_dict[chrom], :])
                
                # diff
                diff_b += self._calc_diff_in_list(beta_old, beta[chrom])
                diff_s += self._calc_diff_in_list(sigma2_old, sigma2[chrom])
            
            diff = diff_b + diff_s
            # niter
            niter += 1
        
        # last update
        # l0, l1 = self._calc_l(yf, ym, h1, h2, beta, sigma2)
        # lld_curr = self._eval_lld(l0, l1, sigma2, n)
        # gamma = self._calc_gamma(l0, l1)
        # lld.append(lld_curr)
        
        out_gamma = {}
        for chrom in chroms:
            out_gamma[chrom] = pd.DataFrame({'prob_z': gamma_list[chrom]})

        return beta, sigma2, out_gamma, lld
        
    def _load_haplotype_from_preloaded(self, path, device):
        return self._to_torch_tensor(np.load(path), device)
        
    def _get_avg_genetic_effect(self, h1, h2, beta):
        return 0.5 * torch.matmul(h1, beta) + 0.5 * torch.matmul(h2, beta)    
        
        
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
    
    def impute(self, df_father, df_mother, df_h1, df_h2, mode, df_covar=None, kwargs={}):
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
        if df_covar is None:
            df_f, df_m, df_1, df_2 = self._extract_and_arrange_rows(
                [df_f, df_m, df_1, df_2],
                individuals_prs
            )
            df_c = None
        else:
            df_f, df_m, df_1, df_2, df_c = self._extract_and_arrange_rows(
                [df_f, df_m, df_1, df_2, df_covar],
                individuals_prs
            )
        
        return impute_method(df_f, df_m, df_1, df_2, covar=df_c, **kwargs)
    
    @staticmethod
    def _get_match_idx(indiv_df, individual_ids):
        target_df = pd.DataFrame({'individual_id': individual_ids})
        indiv_df['idx'] = [ i for i in range(indiv_df.shape[0])]
        target_df = pd.merge(
            target_df,
            indiv_df,
            left_on='individual_id',
            right_on='individual_id',
            how='left'
        )
        return target_df['idx'].values
    
    @staticmethod
    def __is_non_const_col(mat):
        return mat.std(axis=0) != 0
    
    def _get_constant_snp(self, geno_list):
        out = self.__is_non_const_col(geno_list[0])
        if len(geno_list) > 0:
            for i in range(1, len(geno_list)):
                out = np.logical_and(out, self.__is_non_const_col(geno_list[i]))
        return out
    
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
        
        
        

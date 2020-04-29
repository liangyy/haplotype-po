import torch

class BatchLogisticSolver:
    def __init__(self):
        self.info = '''Fit logit(y) ~ x + covariates for every x_i ~ y pair in batch.'''
        self.name = 'BatchLogisticSolver'
        self.stop_criteria = '''Iterate until all x_i ~ y pairs converge'''
    
    def message(self):
        print('-' * 50)
        print('This is {}'.format(self.name))
        print('-' * 50)
        print('* What I do: \n{}'.format(self.info))
        print('* When I stop: \n{}'.format(self.stop_criteria))
    
    def _calc_Mu_S_XSX(self, XSX, X, C, Wcx):
        
        # get k (number of covariates, C)
        k = C.shape[1]
        
        # split Wcx
        Wc = Wcx[:, :k]
        Wx = Wcx[:, k]
        
        # compute w^T x
        C_Wc = torch.einsum('nk,pk->np', C, Wc) 
        X_Wx = torch.einsum('np,p->np', X, Wx)  
        CX_Wcx = C_Wc + X_Wx 
        
        # compute mu and S
        Mu = self._logistic_func(CX_Wcx)
        S = Mu * (1 - Mu)  
        
        # compute X^T S X
        C_S_C = torch.einsum('nk,np,nj->kjp', C, S, C) 
        C_S_X = torch.einsum('nk,np,np->kp', C, S, X)  # .view(k, -1, p)  
        X_S_X = torch.einsum('np,np,np->p', X, S, X)  # .view(1, 1, p)  
        # combine blocks together
        XSX[:k, :k, :] = C_S_C
        XSX[:k, k, :] = C_S_X
        XSX[k, :k, :] = C_S_X
        XSX[k, k, :] = X_S_X
        
        return Mu, S, XSX
    
    def _calc_RHS(self, X, Y, C, Wcx, Mu, XSX):
        
        # get p (number of Xi's)
        p = X.shape[1]
        
        # get RHS
        RES = self._mat_vec_add(-Mu, Y)  
        C_RES = torch.einsum('nk,np->kp', C, RES)  
        X_RES = torch.einsum('np,np->p', X, RES)  
        CX_RES = torch.cat((C_RES, X_RES.view(-1, p)), axis=0)
        # CX_RES[:k, :] = C_RES
        # CX_RES[k, :] = X_RES
        
        XSX_Wcx = torch.einsum('jkp,pk->jp', XSX, Wcx)  
    
        RHS = XSX_Wcx + CX_RES 
        
        return RHS 
    
    @staticmethod
    def _divide_fill_zero(a, b):
        o = torch.div(a, b)
        o[b == 0] = 0
        return o
    
    def _snp_level_convergence_checker(self, w_new, w_old):
        '''
        Return the max and all relative differences between w_new and w_old row-wise.
        '''
        nom = torch.pow(w_new - w_old, 2)
        den = torch.pow(w_new, 2)
        diff = self._divide_fill_zero(torch.sum(nom, axis=1), torch.sum(den, axis=1))
        return diff.max(), diff
    
    def batchIRLS(self, X, y, C, device=None, tol=1e-8, maxiter=100, min_prob=1e-10, use_mask=False):
        '''
        Input
            X: Tensor(n, p)
            y: Tensor(n, 1)
            C: Tensor(n, k)
        Output:
            B_hat: Tensor(p, 1). Effect size estimates where mth entry is for logit(y) ~ X[:, p] + C.
            B_se: Tensor(p, 1). The corresponding SE's of the estimates.
        min_prob: prob < min_prob ~~ 0 and prob > 1 - min_prob ~~ 1 to determine if a node is active or dead
        '''
        if tol > 1 or tol < 0:
            raise ValueError('The tol is relative tolerence so we expect tol in (0, 1).')
            
        # get dimensions
        n = X.shape[0]
        p = X.shape[1]
        k = C.shape[1]
        
        # initialize
        if device is None:
            Wcx = torch.zeros(p, k + 1)
            XSX = torch.Tensor(k + 1, k + 1, p)
            Mu = torch.Tensor(n, p)
            S = torch.Tensor(n, p)
            active_p = torch.zeros(p) == 0
            SE = - torch.ones(k + 1, p)
            diffs = torch.ones(p) 
            SKIP_MASK = torch.ones(p)
        else:
            Wcx = torch.zeros(p, k + 1).to(device)
            XSX = torch.Tensor(k + 1, k + 1, p).to(device)
            Mu = torch.Tensor(n, p).to(device)
            S = torch.Tensor(n, p).to(device)
            active_p = torch.zeros(p).to(device) == 0
            SE = - torch.ones(k + 1, p).to(device)
            diffs = torch.ones(p).to(device) 
            SKIP_MASK = torch.ones(p).to(device) 
        
        max_diff = tol + 1
        # diffs = torch.ones(p) 
        niter = 0
        
        if use_mask is True:
            SKIP_MASK = SKIP_MASK == 0
        else:
            SKIP_MASK = SKIP_MASK == 1
        
        while max_diff > tol and niter < maxiter:
            
            # generate mask
            mask = SKIP_MASK | (diffs > tol)
            
            # take a copy of current Wcx
            Wcx_old = Wcx.clone()
            
            # compute mu, S, and XSX
            Mu[:, mask], _, XSX[:, :, mask] = self._calc_Mu_S_XSX(XSX[:, :, mask], X[:, mask], C, Wcx[mask, :])
            
            # update active status
            active_p[mask] = active_p[mask] & ( (Mu[:, mask] < min_prob).sum(axis=0) == 0 ) & ( (Mu[:, mask] > (1 - min_prob)).sum(axis=0) == 0 )
            mask[mask] = mask[mask] & active_p[mask]
            if mask.sum() == 0:
                break
            
            # get RHS := X^T(S X W + Y - Mu)
            RHS = self._calc_RHS(X[:, mask], y, C, Wcx[mask, :], Mu[:, mask], XSX[:, :, mask]) 

            # solve XSX^{-1} x = RHS and update
            tmp_, _ = torch.solve(
                torch.einsum('kp->pk', RHS).view(mask.sum(), k + 1, -1), 
                torch.einsum('ijk->kij', XSX)[mask, :, :]
            )
            Wcx[mask, :] = tmp_[:, :, 0]
            
            max_diff, diffs[mask] = self._snp_level_convergence_checker(Wcx[mask, :], Wcx_old[mask, :])
            niter += 1
            
        _, _, XSX[:, :, active_p] = self._calc_Mu_S_XSX(XSX[:, :, active_p], X[:, active_p], C, Wcx[active_p, :])
        
        if device is None:
            ONES = torch.eye(k + 1).view(k + 1, k + 1, -1).expand_as(XSX[:, :, active_p])  # Tensor(k + 1, k + 1)
        else:
            ONES = torch.eye(k + 1).to(device).view(k + 1, k + 1, -1).expand_as(XSX[:, :, active_p])   
        
        VAR, _ = torch.solve(
            torch.einsum('ijk->kij', ONES), 
            torch.einsum('ijk->kij', XSX[:, :, active_p])
        )  # Tensor(k + 1, k + 1, p)
        SE[:, active_p] = torch.sqrt(torch.diagonal(VAR, dim1=1, dim2=2)).T 
        
        # return B_hat, B_se
        # if max_diff > tol:
        #     print(f'Warning: not converged! max_diff = {max_diff} > tol = {tol}')
        
        return Wcx.T, SE, diffs < tol
        
    @staticmethod
    def _naive_convergence_checker(w_new, w_old):
        '''
        Return the relative difference between w_new and w_old,
        which is defined as || w_new - w_old ||_2^2 / || w_new ||_2^2.
        '''
        nom = torch.pow(w_new - w_old, 2)
        den = torch.pow(w_new, 2)
        return torch.sum(nom) / torch.sum(den)    
    
    @staticmethod
    def _mat_vec_add(mat, vec):
        '''
        Add vec to each column of mat.
        '''
        nrow = vec.shape[0]
        return mat + vec.view(nrow, 1).expand_as(mat)
        
    @staticmethod
    def _logistic_func(u):
        return 1 / ( 1 + torch.exp(-u) )

# class BatchLogisticSolverWithMask(BatchLogisticSolver):
#     def __init__(self):
#         super().__init__()
#         self.name = 'BatchLogisticSolverWithMask'
#         self.stop_criteria = '''Iterate until x_i ~ y pair converges'''
# 
#     @staticmethod
#     def _divide_fill_zero(a, b):
#         o = torch.div(a, b)
#         o[b == 0] = 0
#         return o
# 
#     def _snp_level_convergence_checker(self, w_new, w_old):
#         '''
#         Return the max and all relative differences between w_new and w_old row-wise.
#         '''
#         nom = torch.pow(w_new - w_old, 2)
#         den = torch.pow(w_new, 2)
#         diff = self._divide_fill_zero(torch.sum(nom, axis=1), torch.sum(den, axis=1))
#         return diff.max(), diff
# 
#     def batchIRLS(self, X, y, C, device=None, tol=1e-8, maxiter=100, min_prob=1e-5):
#         '''
#         Input
#             X: Tensor(n, p)
#             y: Tensor(n, 1)
#             C: Tensor(n, k)
#         Output:
#             B_hat: Tensor(p, 1). Effect size estimates where mth entry is for logit(y) ~ X[:, p] + C.
#             B_se: Tensor(p, 1). The corresponding SE's of the estimates.
#         min_prob: prob < min_prob ~~ 0 and prob > 1 - min_prob ~~ 1 to determine if a node is active or dead
#         '''
#         if tol > 1 or tol < 0:
#             raise ValueError('The tol is relative tolerence so we expect tol in (0, 1).')
# 
#         # get dimensions
#         n = X.shape[0]
#         p = X.shape[1]
#         k = C.shape[1]
# 
#         # initialize
#         if device is None:
#             Wcx = torch.zeros(p, k + 1)
#             XSX = torch.Tensor(k + 1, k + 1, p)
#             Mu = torch.Tensor(n, p)
#             S = torch.Tensor(n, p)
#             active_p = torch.zeros(p) == 0
#             SE = - torch.ones(k + 1, p)
#             diffs = torch.ones(p) 
#         else:
#             Wcx = torch.zeros(p, k + 1).to(device)
#             XSX = torch.Tensor(k + 1, k + 1, p).to(device)
#             Mu = torch.Tensor(n, p).to(device)
#             S = torch.Tensor(n, p).to(device)
#             active_p = torch.zeros(p).to(device) == 1
#             SE = - torch.ones(k + 1, p).to(device)
#             diffs = torch.ones(p).to(device) 
# 
#         max_diff = tol + 1
#         # diffs = torch.ones(p) 
#         niter = 0
# 
#         while max_diff > tol and niter < maxiter:
# 
#             # generate mask
#             mask = diffs > tol 
# 
#             # take a copy of current Wcx
#             Wcx_old = Wcx.clone()
# 
#             # compute mu, S, and XSX
#             Mu[:, mask], _, XSX[:, :, mask] = self._calc_Mu_S_XSX(XSX[:, :, mask], X[:, mask], C, Wcx[mask, :])
# 
#             # update active status
#             active_p[:, mask] = active_p[:, mask] & (Mu[:, mask] > min_prob) & (Mu[:, mask] < 1 - min_prob)
#             mask[:, mask] = mask[:, mask] & active_p[:, mask]
# 
#             # get RHS := X^T(S X W + Y - Mu)
#             RHS = self._calc_RHS(X[:, mask], y, C, Wcx[mask, :], Mu[:, mask], XSX[:, :, mask]) 
# 
#             # solve XSX^{-1} x = RHS and update
#             tmp_, _ = torch.solve(
#                 torch.einsum('kp->pk', RHS).view(mask.sum(), k + 1, -1), 
#                 torch.einsum('ijk->kij', XSX)[mask, :, :]
#             )
#             Wcx[mask, :] = tmp_[:, :, 0]
# 
#             max_diff, diffs[mask] = self._snp_level_convergence_checker(Wcx[mask, :], Wcx_old[mask, :])
#             niter += 1
# 
#         _, _, XSX[:, :, active_p] = self._calc_Mu_S_XSX(XSX[:, :, active_p], X[:, active_p], C, Wcx[active_p, :])
# 
#         if device is None:
#             ONES = torch.eye(k + 1).view(k + 1, k + 1, -1).expand_as(XSX[:, :, active_p])  # Tensor(k + 1, k + 1)
#         else:
#             ONES = torch.eye(k + 1).to(device).view(k + 1, k + 1, -1).expand_as(XSX[:, :, active_p])   
# 
#         VAR, _ = torch.solve(
#             torch.einsum('ijk->kij', ONES), 
#             torch.einsum('ijk->kij', XSX[:, :, active_p])
#         )  # Tensor(k + 1, k + 1, p)
#         SE[:, active_p] = torch.sqrt(torch.diagonal(VAR, dim1=1, dim2=2)).T 
# 
#         # return B_hat, B_se
#         # if max_diff > tol:
#         #     print(f'Warning: not converged! max_diff = {max_diff} > tol = {tol}')
# 
# 
#         return Wcx.T, SE, diffs < tol
# 

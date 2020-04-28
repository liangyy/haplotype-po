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
    
    def batchIRLS(self, X, y, C, device=None, tol=1e-8, maxiter=100):
        '''
        Input
            X: Tensor(n, p)
            y: Tensor(n, 1)
            C: Tensor(n, k)
        Output:
            B_hat: Tensor(p, 1). Effect size estimates where mth entry is for logit(y) ~ X[:, p] + C.
            B_se: Tensor(p, 1). The corresponding SE's of the estimates.
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
            XSX = torch.Tensor(k + 1, k + 1, pe)
            diffs = torch.ones(p) 
        else:
            Wcx = torch.zeros(p, k + 1).to(device)
            XSX = torch.Tensor(k + 1, k + 1, p).to(device)
            diffs = torch.ones(p).to(device) 
        
        diff = tol + 1
        niter = 0
        
        while diff > tol and niter < maxiter:
            
            # take a copy of current Wcx
            Wcx_old = Wcx.clone()
            
            # compute mu, S, and X^T S X
            Mu, S, XSX = self._calc_Mu_S_XSX(XSX, X, C, Wcx)
            
            # get RHS := X^T(S X W + Y - Mu)
            RHS = self._calc_RHS(X, y, C, Wcx, Mu, XSX)

            # solve XSX^{-1} x = RHS and update
            Wcx, _ = torch.solve(
                torch.einsum('kp->pk', RHS).view(p, k + 1, -1), 
                torch.einsum('ijk->kij', XSX)
            )
            Wcx = Wcx[:, :, 0]
            
            diff = self._naive_convergence_checker(Wcx, Wcx_old)
            niter += 1
        
        # after iteration
        _, _, XSX = self._calc_Mu_S_XSX(XSX, X, C, Wcx)
        
        if device is None:
            ONES = torch.eye(k + 1).view(k + 1, k + 1, -1).expand_as(XSX)  # Tensor(k + 1, k + 1)
        else:
            ONES = torch.eye(k + 1).to(device).view(k + 1, k + 1, -1).expand_as(XSX)  
            
        VAR, _ = torch.solve(
            torch.einsum('ijk->kij', ONES), 
            torch.einsum('ijk->kij', XSX)
        )  # Tensor(k + 1, k + 1, p)
        
        # return B_hat, B_se
        if diff > tol:
            print(f'Warning: not converged! diff = {diff} > tol = {tol}')
        
        return Wcx.T, torch.sqrt(torch.diagonal(VAR, dim1=1, dim2=2)).T 
        
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

class BatchLogisticSolverWithMask(BatchLogisticSolver):
    def __init__(self):
        super().__init__()
        self.name = 'BatchLogisticSolverWithMask'
        self.stop_criteria = '''Iterate until x_i ~ y pair converges'''
    
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
    
    def batchIRLS(self, X, y, C, device=None, tol=1e-8, maxiter=100):
        '''
        Input
            X: Tensor(n, p)
            y: Tensor(n, 1)
            C: Tensor(n, k)
        Output:
            B_hat: Tensor(p, 1). Effect size estimates where mth entry is for logit(y) ~ X[:, p] + C.
            B_se: Tensor(p, 1). The corresponding SE's of the estimates.
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
            diffs = torch.ones(p) 
        else:
            Wcx = torch.zeros(p, k + 1).to(device)
            XSX = torch.Tensor(k + 1, k + 1, p).to(device)
            diffs = torch.ones(p).to(device) 
        
        max_diff = tol + 1
        # diffs = torch.ones(p) 
        niter = 0
        
        while max_diff > tol and niter < maxiter:
            
            # generate mask
            mask = diffs > tol
            
            # take a copy of current Wcx
            Wcx_old = Wcx.clone()
            
            # compute mu, S, and XSX
            Mu, S, XSX[:, :, mask] = self._calc_Mu_S_XSX(XSX[:, :, mask], X[:, mask], C, Wcx[mask, :])
            
            # get RHS := X^T(S X W + Y - Mu)
            RHS = self._calc_RHS(X[:, mask], y, C, Wcx[mask, :], Mu, XSX[:, :, mask]) 

            # solve XSX^{-1} x = RHS and update
            tmp_, _ = torch.solve(
                torch.einsum('kp->pk', RHS).view(mask.sum(), k + 1, -1), 
                torch.einsum('ijk->kij', XSX)[mask, :, :]
            )
            Wcx[mask, :] = tmp_[:, :, 0]
            
            max_diff, diffs = self._snp_level_convergence_checker(Wcx, Wcx_old)
            niter += 1
            
        _, _, XSX = self._calc_Mu_S_XSX(XSX, X, C, Wcx)
        
        if device is None:
            ONES = torch.eye(k + 1).view(k + 1, k + 1, -1).expand_as(XSX)  # Tensor(k + 1, k + 1)
        else:
            ONES = torch.eye(k + 1).to(device).view(k + 1, k + 1, -1).expand_as(XSX)   
        
        VAR, _ = torch.solve(
            torch.einsum('ijk->kij', ONES), 
            torch.einsum('ijk->kij', XSX)
        )  # Tensor(k + 1, k + 1, p)
        
        # return B_hat, B_se
        if max_diff > tol:
            print(f'Warning: not converged! max_diff = {max_diff} > tol = {tol}')
        return Wcx.T, torch.sqrt(torch.diagonal(VAR, dim1=1, dim2=2)).T 
        

import xarray as xr
from .log_progress import log_progress
import numpy as np

def get_neibours(xi,eta,r,grid):
    neibours = []
    r=int(r/2)
    for n_xi in np.arange(0,r+1):
        for n_eta in np.arange(0,r+1):
            if (n_xi**2 + n_eta**2) > r**2:
                break
            else:
                if (xi-n_xi >= 0) and (eta-n_eta >= 0):
                    neibours.append([xi-n_xi,eta-n_eta])
                elif (xi+n_xi <= grid.xi_rho[-1]) and (eta+n_eta <= grid.eta_rho[-1]):
                    neibours.append([xi+n_xi,eta+n_eta])
    return neibours

def mask_radius(mask_da,r):
    mask_new = mask_da.copy()
    for xi in log_progress(mask_da.xi_rho.values,name='xi'):
        for eta in mask_da.eta_rho.values:
            if mask_da.isel(eta_rho=eta,xi_rho=xi) != 0:
                neibours = get_neibours(xi,eta,r,mask_da)
                for neib in neibours:
                    if mask_da.isel(xi_rho=neib[0],eta_rho=neib[1]).values == 0:
                        mask_new[eta,xi]=0
                        break
    return mask_new
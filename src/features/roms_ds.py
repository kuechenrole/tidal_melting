#read in raw data as xr.dataset
import xarray as xr
import numpy as np
from .calc_z import calc_z

def make_4D_mask(ds):

    
    mask_4d = np.tile(ds.mask_rho,(ds.ocean_time.size,ds.s_rho.size,1,1))
    ds['mask_4d'] = xr.DataArray(mask_4d,dims=['ocean_time','s_rho','eta_rho','xi_rho'])
    ds.mask_4d.attrs = ds.mask_rho.attrs
    
    return ds



def make_3D_XiEta(ds):
    
    xi_3d = np.tile(ds.xi_rho,(ds.s_rho.size,ds.eta_rho.size,1))
    eta_3d = np.swapaxes(np.tile(ds.eta_rho,(ds.s_rho.size,ds.xi_rho.size,1)),1,2)
    
    xi_3d_da = xr.DataArray(xi_3d,dims=['s_rho','eta_rho','xi_rho'])
    eta_3d_da = xr.DataArray(eta_3d,dims=['s_rho','eta_rho','xi_rho'])
    
    ds = ds.assign_coords(xi_3d=xi_3d_da)
    ds = ds.assign_coords(eta_3d=eta_3d_da)
    
    ds['xi_3d'] = ds.xi_3d.where(ds.mask_rho == 1)
    ds['eta_3d'] = ds.eta_3d.where(ds.mask_rho ==1)
    
    ds.xi_3d.attrs = ds.xi_rho.attrs
    ds.eta_3d.attrs = ds.eta_rho.attrs
    
    return ds



def make_4D_depth(ds):
    
    depths = np.empty((ds.ocean_time.size,ds.s_rho.size,ds.eta_rho.size,ds.xi_rho.size))
    
    for tstep in np.arange(ds.ocean_time.size):

        h = ds.h.values
        zice = ds.zice.values
        theta_s = ds.theta_s.values
        theta_b = ds.theta_b.values
        hc = ds.hc.values
        N = ds.s_rho.size
        zeta = ds.zeta[tstep].values
        Vstretching = ds.Vstretching.values
        
        depths[tstep],s,C = calc_z(h,zice,theta_s,theta_b,hc,N,zeta,Vstretching)
        
    ds = ds.assign_coords(depth = xr.DataArray(depths,dims=['ocean_time','s_rho','eta_rho','xi_rho']))
    
    ds['depth'] = ds.depth.where(ds.mask_rho == 1)
    
    return ds



def make_roms_ds(file_paths):
    '''Takes a roms history or averages file (wildcards are possible) and returns a Xarray dataset including 4D mask, 3D grid coordinates and 4D depths'''
    
    print('set up multifile dataset')
    ds_tmp = xr.open_mfdataset(file_paths,data_vars='minimal')
    
    print('set up 4D mask and add as variable to dataset')
    ds_tmp = make_4D_mask(ds_tmp)
    
    print('set up 3D xi and eta arrays, fill with NaNs where invalid and apply as coordinates')
    ds_tmp = make_3D_XiEta(ds_tmp)
    
    print('calculate 4D depth array, fill with NaNs where invalid and apply as coordinate')
    ds = make_4D_depth(ds_tmp)
    
    return ds
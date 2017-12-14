#read in raw data as xr.dataset
import xarray as xr
import numpy as np
from .calc_z import calc_z

def make_roms_ds(file_paths):
    
    print('set up multifile dataset')
    ds_tmp = xr.open_mfdataset(file_paths)
    
    
    print('set up 4D mask and add as variable to dataset')
    mask_4d = np.swapaxes(np.tile(ds_tmp.mask_rho,(ds_tmp.s_rho.size,1,1,1)),0,1)
    mask_4d_da = xr.DataArray(mask_4d,dims=['ocean_time','s_rho','eta_rho','xi_rho'])
    ds_tmp['mask_4d'] = mask_4d_da
    ds_tmp.mask_4d.attrs = ds_tmp.mask_rho.attrs
    
    
    print('set up 3D xi and eta arrays, mask them and apply as coordinates')
    xi_3d = np.tile(ds_tmp.xi_rho,(ds_tmp.s_rho.size,ds_tmp.eta_rho.size,1))
    eta_3d = np.swapaxes(np.tile(ds_tmp.eta_rho,(ds_tmp.s_rho.size,ds_tmp.xi_rho.size,1)),1,2)
    
    xi_3d_ma = np.ma.masked_where(mask_4d[0]==0,xi_3d)
    eta_3d_ma = np.ma.masked_where(mask_4d[0]==0,eta_3d)
    
    xi_3d_da = xr.DataArray(xi_3d_ma,dims=['s_rho','eta_rho','xi_rho'])
    eta_3d_da = xr.DataArray(eta_3d_ma,dims=['s_rho','eta_rho','xi_rho'])
    
    ds_tmp = ds_tmp.assign_coords(xi_3d=xi_3d_da)
    ds_tmp = ds_tmp.assign_coords(eta_3d=eta_3d_da)
    
    ds_tmp.xi_3d.attrs = ds_tmp.xi_rho.attrs
    ds_tmp.eta_3d.attrs = ds_tmp.eta_rho.attrs
    
    
    print('calculate 4D depth array')
    depths = np.empty((ds_tmp.ocean_time.size,ds_tmp.s_rho.size,ds_tmp.eta_rho.size,ds_tmp.xi_rho.size))
    
    for tstep in np.arange(ds_tmp.ocean_time.size):

        h = ds_tmp.h[tstep].values
        zice = ds_tmp.zice[tstep].values
        theta_s = ds_tmp.theta_s[tstep].values
        theta_b = ds_tmp.theta_b[tstep].values
        hc = ds_tmp.hc[tstep].values
        N = ds_tmp.s_rho.size
        zeta = ds_tmp.zeta[tstep].values
        Vstretching = ds_tmp.Vstretching[tstep].values
        depths[tstep],s,C = calc_z(h,zice,theta_s,theta_b,hc,N,zeta,Vstretching)
        
    print('apply mask to depths')
    depths_ma = np.ma.masked_where(ds_tmp.mask_4d==0,depths)
        
    print('assign depth as new coordinate to the data set')
    depth_da = xr.DataArray(depths_ma,dims=['ocean_time','s_rho','eta_rho','xi_rho'])
    ds = ds_tmp.assign_coords(depth=depth_da)
    
    return ds
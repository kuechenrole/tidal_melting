#read in raw data as xr.dataset
import xarray as xr
import numpy as np
from .calc_z import calc_z
from .cartesian_grid_3d import cartesian_grid_3d
from .rotate_vector_roms import rotate_vector_roms

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
    
    #ds['xi_3d'] = ds.xi_3d.where(ds.mask_rho == 1)
    #ds['eta_3d'] = ds.eta_3d.where(ds.mask_rho ==1)
    
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
    
    #ds['depth'] = ds.depth.where(ds.mask_rho == 1)
    
    return ds

def make_cartesian_grid_3D(ds):
    
    lon_u = ds.lon_u.values
    lat_u = ds.lat_u.values
    lon_v = ds.lon_v.values
    lat_v = ds.lat_v.values
    h = ds.h.values
    zice = ds.zice.values
    theta_s = ds.theta_s.values
    theta_b = ds.theta_b.values
    hc = ds.hc.values
    N = ds.s_rho.size
    
    z = np.empty((ds.ocean_time.size,ds.s_rho.size,ds.eta_rho.size,ds.xi_rho.size))
    dz = np.empty(np.shape(z))
    
    for tstep in np.arange(ds.ocean_time.size):
        
        zeta = ds.zeta[tstep].values
        
        dx,dy,dz[tstep],z[tstep] = cartesian_grid_3d(lon_u, lat_u, lon_v, lat_v, h, zice, theta_s, theta_b, hc, N, zeta)
        
    ds['dx'] = xr.DataArray(dx,dims=['s_rho','eta_rho','xi_rho'])
    ds['dx'] = ds.dx.where(ds.mask_rho == 1)
    
    ds['dy'] = xr.DataArray(dy,dims=['s_rho','eta_rho','xi_rho'])
    ds['dy'] = ds.dy.where(ds.mask_rho == 1)
    
    ds['dz'] = xr.DataArray(dz,dims=['ocean_time','s_rho','eta_rho','xi_rho'])
    ds['dz'] = ds.dz.where(ds.mask_rho == 1)
    
    ds['z'] = xr.DataArray(z,dims=['ocean_time','s_rho','eta_rho','xi_rho'])
    ds['z'] = ds.z.where(ds.mask_rho == 1)
    
    dV = ds.dx * ds.dy * ds.dz
    ds['dV'] = dV
    
    return ds

def make_uvbar_lonlat(ds):
    
    angle = ds.angle.values
    
    ubar_lonlat = np.empty((ds.ocean_time.size,ds.eta_rho.size,ds.xi_rho.size))
    vbar_lonlat = np.empty((ds.ocean_time.size,ds.eta_rho.size,ds.xi_rho.size))
    
    for tstep in np.arange(ds.ocean_time.size):
        ubar=ds.ubar[tstep].values
        vbar=ds.vbar[tstep].values
        
        ubar_lonlat[tstep],vbar_lonlat[tstep] = rotate_vector_roms(ubar,vbar,angle)
        
    ds['ubar_lonlat'] = xr.DataArray(ubar_lonlat,dims=['ocean_time','eta_rho','xi_rho'])
    ds['ubar_lonlat'] = ds.ubar_lonlat.where(ds.mask_rho == 1) 
    ds.ubar_lonlat.attrs = ds.u_bar.attrs
    
    ds['vbar_lonlat'] = xr.DataArray(vbar_lonlat,dims=['ocean_time','eta_rho','xi_rho'])
    ds['vbar_lonlat'] = ds.vbar_lonlat.where(ds.mask_rho == 1)
    ds.vbar_lonlat.attrs = ds.vbar.attrs
    
    return ds

def make_uv_lonlat(ds):
    
    angle = ds.angle.values
    
    u_lonlat = np.empty((ds.ocean_time.size,ds.s_rho.size,ds.eta_rho.size,ds.xi_rho.size))
    v_lonlat = np.empty((ds.ocean_time.size,ds.s_rho.size,ds.eta_rho.size,ds.xi_rho.size))
    
    for tstep in np.arange(ds.ocean_time.size):
        for level in np.arange(ds.s_rho.size):
            u=ds.u[tstep,level].values
            v=ds.v[tstep,level].values
        
            u_lonlat[tstep,level],v_lonlat[tstep,level] = rotate_vector_roms(u,v,angle)
        
    ds['u_lonlat'] = xr.DataArray(u_lonlat,dims=['ocean_time','s_rho','eta_rho','xi_rho'])
    ds['u_lonlat'] = ds.u_lonlat.where(ds.mask_rho == 1) 
    ds.u_lonlat.attrs = ds.u.attrs
    
    ds['v_lonlat'] = xr.DataArray(v_lonlat,dims=['ocean_time','s_rho','eta_rho','xi_rho'])
    ds['v_lonlat'] = ds.v_lonlat.where(ds.mask_rho == 1)
    ds.v_lonlat.attrs = ds.v.attrs
    
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
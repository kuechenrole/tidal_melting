import xarray as xr
import numpy as np
src_dir = os.path.join(os.environ.get('projdir'),'src')
sys.path.append(src_dir)
from features.rotate_vector_roms import rotate_vector_roms
from features.calc_z import calc_z
import gsw

def make_cartesian_grid_3D(grd,ds):

    lon_u = grd.lon_u.values
    lat_u = grd.lat_u.values
    lon_v = grd.lon_v.values
    lat_v = grd.lat_v.values
    h = grd.h.values
    zice = grd.zice.values
    theta_s = ds.theta_s.values
    theta_b = ds.theta_b.values
    hc = ds.hc.values
    N = ds.s_rho.size

    zeta = ds.zeta.values

    dx,dy,dz,z = cartesian_grid_3d(lon_u, lat_u, lon_v, lat_v, h, zice, theta_s, theta_b, hc, N, zeta)

    ds['dx'] = xr.DataArray(dx,dims=['s_rho','eta_rho','xi_rho'])
    ds['dx'] = ds.dx.where(ds.mask_rho == 1)

    ds['dy'] = xr.DataArray(dy,dims=['s_rho','eta_rho','xi_rho'])
    ds['dy'] = ds.dy.where(ds.mask_rho == 1)

    ds['dz'] = xr.DataArray(dz,dims=['s_rho','eta_rho','xi_rho'])
    ds['dz'] = ds.dz.where(ds.mask_rho == 1)

    ds['z'] = xr.DataArray(z,dims=['s_rho','eta_rho','xi_rho'])
    ds['z'] = ds.z.where(ds.mask_rho == 1)

    dV = ds.dx * ds.dy * ds.dz
    ds['dV'] = dV

    return ds

def make_uv_lonlat(grd,ds):

    angle = grd.angle.values

    u_lonlat = np.empty((ds.s_rho.size,ds.eta_rho.size,ds.xi_rho.size))
    v_lonlat = np.empty((ds.s_rho.size,ds.eta_rho.size,ds.xi_rho.size))


    for level in np.arange(ds.s_rho.size):
        u=ds.u[level].values
        v=ds.v[level].values

        u_lonlat[level],v_lonlat[level] = rotate_vector_roms(u,v,angle)

    ds['u_lonlat'] = xr.DataArray(u_lonlat,dims=['s_rho','eta_rho','xi_rho'])
    ds['u_lonlat'] = ds.u_lonlat.where(ds.mask_rho == 1)
    ds.u_lonlat.attrs = ds.u.attrs

    ds['v_lonlat'] = xr.DataArray(v_lonlat,dims=['s_rho','eta_rho','xi_rho'])
    ds['v_lonlat'] = ds.v_lonlat.where(ds.mask_rho == 1)
    ds.v_lonlat.attrs = ds.v.attrs

    return ds

def make_depth(grd,ds):

    h = grd.h.values
    zice = grd.zice.values
    theta_s = ds.theta_s.values
    theta_b = ds.theta_b.values
    hc = ds.hc.values
    N = ds.s_rho.size
    zeta = ds.zeta.values
    Vstretching = ds.Vstretching.values

    depths,s,C = calc_z(h,zice,theta_s,theta_b,hc,N,zeta,Vstretching)

    ds = ds.assign_coords(depth = xr.DataArray(depths,dims=['s_rho','eta_rho','xi_rho']))

    return ds

def make_density(grd,ds):
    
    p = gsw.conversions.p_from_z(ds.depth,grd.lat_rho)
    SA = gsw.conversions.SA_from_SP(ds.salt,p,grd.lon_rho,grd.lat_rho)
    CT = gsw.conversions.CT_from_pt(SA,ds.temp)
    
    ds['rho']=(('s_rho','eta_rho','xi_rho'),gsw.density.rho(SA,CT,p))
    
    return ds
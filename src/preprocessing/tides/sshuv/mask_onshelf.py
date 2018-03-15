from scipy.ndimage.filters import gaussian_filter, uniform_filter
import os
import xarray as xr

out_path = os.path.join(os.environ.get('prodir'),'waom5_tds.nc')
tds_path = os.path.join(os.environ.get('intdir'),'waom5_tds.nc')
grid_path =  os.path.join(os.environ.get('prodir'),'waom5_grd.nc')


tds_ds = xr.open_dataset(tds_path)

grid_ds = xr.open_dataset(grid_path)

mask_tmp = grid_ds.mask_rho.where((grid_ds.zice == 0),0.0)
mask_tmp = mask_tmp.where(grid_ds.h > 1000,0.0)
mask_tmp.values = gaussian_filter(mask_tmp,20)

tds_ds['tide_Pamp']=tds_ds.tide_Pamp*mask_tmp
tds_ds.to_netcdf(out_path,'w')

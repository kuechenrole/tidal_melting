import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from shutil import copyfile

int_grid_path = os.path.join(os.environ.get('intdir'),'waom5_grd.nc')
pro_grid_path = os.path.join(os.environ.get('prodir'),'waom5_grd.nc')

copyfile(int_grid_path,pro_grid_path)

grid = xr.open_dataset(int_grid_path)

grid['zice'] = grid.zice.where(grid.lat_rho<-60.0,0.0)
grid.zice[736,394]=0.0

grid.zice.to_dataset().to_netcdf(pro_grid_path,'a')


wct = grid.h + grid.zice
grid['h'] = grid.h.where((wct.values>=30.0) | (grid.zice.values==0.0) ,-grid.zice+30.0)
wct_new = grid.h + grid.zice
wct_new[943,444]

grid.h.to_dataset().to_netcdf(pro_grid_path,'a')

#grid.to_netcdf(pro_grid_path)
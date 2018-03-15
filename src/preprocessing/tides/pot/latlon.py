import numpy as np
import netCDF4
import sys

# Open the grid file and read lon, lat
ncfile = '/home/ubuntu/bigStick/waom10Grids/waom10_grd_medium.nc'
nc = netCDF4.Dataset(ncfile, 'r', format='NETCDF3_CLASSIC')
lon_rho = nc.variables['lon_rho'][:]
lat_rho = nc.variables['lat_rho'][:]
nc.close()

# Print out the requested information
dims = lon_rho.shape
print(dims[0]*dims[1])

for j in range(dims[0]):
    for i in range(dims[1]):
        if lon_rho[j,i]<0:
            lon_rho[j,i]=lon_rho[j,i]+360.0
        print(str(lon_rho[j,i]) + "   " + str(lat_rho[j,i]))


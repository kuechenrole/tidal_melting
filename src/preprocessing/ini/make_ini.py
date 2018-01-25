import xarray as xr
from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import NearestNDInterpolator
import sys
import os
sys.path.append(os.path.join(os.environ['projdir'],'src'))
from features.calc_z import *


# Main routine
def roms_ini_ecco(grid_file, theta_file, salt_file, output_file, Tcline, theta_s, theta_b, hc, N, nbdry_ecco,tempbounds,saltbounds,vostock):

    # Read ECCO2 data and grid
    print('Reading ECCO2 data')
    theta_ds = xr.open_dataset(theta_file)
    theta_raw = transpose(theta_ds.THETA.values[0,:,0:nbdry_ecco,:])
    lon_ecco_raw = theta_ds.LONGITUDE_T.values[:]
    lat_ecco = theta_ds.LATITUDE_T.values[0:nbdry_ecco]
    depth_ecco_raw = theta_ds.DEPTH_T.values[:]

    salt_ds = xr.open_dataset(salt_file)
    salt_raw = transpose(salt_ds.SALT.values[0,:,0:nbdry_ecco,:])
    # The ECCO2 longitude axis doesn't wrap around; there is a gap between
    # almost-180W and almost-180E, and the ROMS grid has points in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_ecco = zeros(size(lon_ecco_raw)+2)
    lon_ecco[0] = lon_ecco_raw[-1]-360
    lon_ecco[1:-1] = lon_ecco_raw
    lon_ecco[-1] = lon_ecco_raw[0]+360

    # The shallowest ECCO2 depth value is 5 m, but ROMS needs 0 m. So add the
    # index depth = 0 m to the beginning. Later we will just copy the 5 m values
    # for theta and salt into this index. Similarly, add the index depth = 6000 m
    # to the end.
    depth_ecco = zeros(size(depth_ecco_raw)+2)
    depth_ecco[0] = 0.0
    depth_ecco[1:-1] = depth_ecco_raw
    depth_ecco[-1] = 6000.0

    # Copy the theta and salt values to the new longitude and depth indices,
    # making sure to preserve the mask.
    theta = ma.array(zeros((size(lon_ecco), size(lat_ecco), size(depth_ecco))))
    theta[1:-1,:,1:-1] = ma.copy(theta_raw)
    theta[0,:,1:-1] = ma.copy(theta_raw[-1,:,:])
    theta[-1,:,1:-1] = ma.copy(theta_raw[0,:,:])
    theta[:,:,0] = ma.copy(theta[:,:,1])
    theta[:,:,-1] = ma.copy(theta[:,:,-2])
    salt = ma.array(zeros((size(lon_ecco), size(lat_ecco), size(depth_ecco))))
    salt[1:-1,:,1:-1] = ma.copy(salt_raw)
    salt[0,:,1:-1] = ma.copy(salt_raw[-1,:,:])
    salt[-1,:,1:-1] = ma.copy(salt_raw[0,:,:])
    salt[:,:,0] = ma.copy(salt[:,:,1])
    salt[:,:,-1] = ma.copy(salt[:,:,-2])

    # Read ROMS grid
    print('Reading ROMS grid')
    grid_fid = Dataset(grid_file, 'r')
    lon_roms = grid_fid.variables['lon_rho'][:,:]
    #ecco goes from 0 to 360 roms from -180 to 180
    for i in range(len(lon_roms)):
        for j in range(len(lon_roms[i])):
            if lon_roms[i,j] < 0:
                lon_roms[i,j] = lon_roms[i,j] + 360
    lat_roms = grid_fid.variables['lat_rho'][:,:]
    h = grid_fid.variables['h'][:,:]
    zice = grid_fid.variables['zice'][:,:]
    mask_rho = grid_fid.variables['mask_rho'][:,:]
    #mask_zice = grid_fid.variables['mask_zice'][:,:]
    grid_fid.close()
    mask_zice = where(zice < 0.0,1,0)*mask_rho
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)
    # Mask h and zice with zeros
    h = h*mask_rho
    zice = zice*mask_zice

    # Get a 3D array of ROMS z-coordinates, as well as 1D arrays of s-coordinates
    # and stretching curves
    z_roms_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N)
    # Copy the latitude and longitude values into 3D arrays of the same shape
    lon_roms_3d = tile(lon_roms, (N,1,1))
    lat_roms_3d = tile(lat_roms, (N,1,1))

    # Regridding happens here...
    print('Interpolating temperature')
    temp = interp_ecco2roms(theta, lon_ecco, lat_ecco, depth_ecco, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice,-0.5,vostock)
    print('Interpolating salinity')
    salt = interp_ecco2roms(salt, lon_ecco, lat_ecco, depth_ecco, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice,34.5,vostock)

    # setting bounds for temperture and salinity
    tempLow,tempUp=tempbounds
    saltLow,saltUp=saltbounds
    temp[temp < tempLow] = tempLow
    temp[temp > tempUp] = tempUp
    salt[salt < saltLow] = saltLow
    salt[salt > saltUp] = saltUp 

    # Set initial velocities and sea surface height to zero
    u = zeros((N, num_lat, num_lon-1))
    v = zeros((N, num_lat-1, num_lon))
    ubar = zeros((num_lat, num_lon-1))
    vbar = zeros((num_lat-1, num_lon))
    zeta = zeros((num_lat, num_lon))

    print('Writing to NetCDF file')
    out_fid = Dataset(output_file, 'w')
    # Define dimensions
    out_fid.createDimension('xi_u', num_lon-1)
    out_fid.createDimension('xi_v', num_lon)
    out_fid.createDimension('xi_rho', num_lon)
    out_fid.createDimension('eta_u', num_lat)
    out_fid.createDimension('eta_v', num_lat-1)
    out_fid.createDimension('eta_rho', num_lat)
    out_fid.createDimension('s_rho', N)
    out_fid.createDimension('ocean_time', None)
    out_fid.createDimension('one', 1);
    # Define variables and assign values
    out_fid.createVariable('tstart', 'f8', ('one'))
    out_fid.variables['tstart'].long_name = 'start processing day'
    out_fid.variables['tstart'].units = 'day'
    out_fid.variables['tstart'][:] = 0.0
    out_fid.createVariable('tend', 'f8', ('one'))
    out_fid.variables['tend'].long_name = 'end processing day'
    out_fid.variables['tend'].units = 'day'
    out_fid.variables['tend'][:] = 0.0
    out_fid.createVariable('theta_s', 'f8', ('one'))
    out_fid.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
    out_fid.variables['theta_s'][:] = theta_s
    out_fid.createVariable('theta_b', 'f8', ('one'))
    out_fid.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
    out_fid.variables['theta_b'].units = 'nondimensional'
    out_fid.variables['theta_b'][:] = theta_b
    out_fid.createVariable('Tcline', 'f8', ('one'))
    out_fid.variables['Tcline'].long_name = 'S-coordinate surface/bottom layer width'
    out_fid.variables['Tcline'].units = 'meter'
    out_fid.variables['Tcline'][:] = Tcline
    out_fid.createVariable('hc', 'f8', ('one'))
    out_fid.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
    out_fid.variables['hc'].units = 'meter'
    out_fid.variables['hc'][:] = hc
    out_fid.createVariable('Cs_r', 'f8', ('s_rho'))
    out_fid.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    out_fid.variables['Cs_r'].units = 'nondimensional'
    out_fid.variables['Cs_r'].valid_min = -1.0
    out_fid.variables['Cs_r'].valid_max = 0.0
    out_fid.variables['Cs_r'][:] = Cs_r
    out_fid.createVariable('ocean_time', 'f8', ('ocean_time'))
    out_fid.variables['ocean_time'].long_name = 'time since initialization'
    out_fid.variables['ocean_time'].units = 'seconds'
    out_fid.variables['ocean_time'][0] = 0.0
    out_fid.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'))
    out_fid.variables['u'].long_name = 'u-momentum component'
    out_fid.variables['u'].units = 'meter second-1'
    out_fid.variables['u'][0,:,:,:] = u
    out_fid.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'))
    out_fid.variables['v'].long_name = 'v-momentum component'
    out_fid.variables['v'].units = 'meter second-1'
    out_fid.variables['v'][0,:,:,:] = v
    out_fid.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'))
    out_fid.variables['ubar'].long_name = 'vertically integrated u-momentum component'
    out_fid.variables['ubar'].units = 'meter second-1'
    out_fid.variables['ubar'][0,:,:] = ubar
    out_fid.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'))
    out_fid.variables['vbar'].long_name = 'vertically integrated v-momentum component'
    out_fid.variables['vbar'].units = 'meter second-1'
    out_fid.variables['vbar'][0,:,:] = vbar
    out_fid.createVariable('zeta', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'))
    out_fid.variables['zeta'].long_name = 'free-surface'
    out_fid.variables['zeta'].units = 'meter'
    out_fid.variables['zeta'][0,:,:] = zeta
    out_fid.createVariable('temp', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    out_fid.variables['temp'].long_name = 'potential temperature'
    out_fid.variables['temp'].units = 'Celsius'
    out_fid.variables['temp'][0,:,:,:] = temp
    out_fid.createVariable('salt', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    out_fid.variables['salt'].long_name = 'salinity'
    out_fid.variables['salt'].units = 'PSU'
    out_fid.variables['salt'][0,:,:,:] = salt
    out_fid.createVariable('sc_r', 'f8', ('s_rho'))
    out_fid.variables['sc_r'].long_name = 'S-coordinate at rho-points'
    out_fid.variables['sc_r'].units = 'nondimensional'
    out_fid.variables['sc_r'].valid_min = -1.0
    out_fid.variables['sc_r'].valid_max = 0.0
    out_fid.variables['sc_r'][:] = sc_r
    out_fid.close()

    print('Finished')

def interp_ecco2roms (A, lon_ecco, lat_ecco, depth_ecco, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice,fill,vostock):

    # Calculate N based on size of ROMS grid
    N = size(lon_roms_3d, 0)

    eta_roms,xi_roms=mgrid[0:size(lon_roms_3d,1),0:size(lon_roms_3d,2)]

    # Unmask A and fill with NaNs
    A_unmask = A.data
    A_unmask[A.mask]=NaN

    # Build a function for linear interpolation of A; set flag to fill
    # out-of-bounds values with NaN
    interp_function = RegularGridInterpolator((lon_ecco, lat_ecco, depth_ecco), A_unmask, bounds_error=False, fill_value=NaN)
    B = zeros(shape(lon_roms_3d))

    # Interpolate each z-level individually - 3D vectorisation uses too
    # much memory!
    for k in range(N):
        print('...vertical level ', str(k+1), ' of ', str(N))
        # Pass positive values for ROMS depth
        B[k] = interp_function((lon_roms_3d[k,:,:], lat_roms_3d[k,:,:], -z_roms_3d[k,:,:]),method='nearest')
        # Fill NaNs with constant value
        #index = isnan(tmp)
        #tmp[index] = fill
        # Fill land mask with constant value
        #tmp[mask_rho==0] = fill
        # Save this depth level
        #B[k,:,:] = tmp
    print('Extrapolating under ice shelves')   
    
    for k in range(N):
            
        print('processing depth level: ',k)
        data = B[k]

        valid_mask = ~isnan(data)
        coords = array(nonzero(valid_mask)).T
        values = data[valid_mask]

        it = NearestNDInterpolator(coords,values)

        filled = it(list(ndindex(data.shape))).reshape(data.shape)

        # Fill in known values on high res grid
        B[k] = filled
    

    return B

# set user parameter and call main routine

grd_file = os.path.join(os.environ['projdir'],'data','preprocessing','processed','waom10_grd.nc')
out_file = os.path.join(os.environ['projdir'],'data','preprocessing','processed','waom10_ini.nc')

print('making ini file for grid file: ' + grd_file)

# Path to ECCO2 files for temperature and salinity 1th January 1995
salt_file = os.path.join(os.environ['projdir'],'data','preprocessing','external','ecco2','SALT.nc','SALT.1440x720x50.20070101.nc')
theta_file = os.path.join(os.environ['projdir'],'data','preprocessing','external','ecco2','THETA.nc','THETA.1440x720x50.20070101.nc')

# Grid parameters; check grid_file and *.in file to make sure these are correct
Tcline = 50
theta_s = 4.0
theta_b = 0.9
hc = 50
N = 31

# Northernmost index of ECCO2 grid to read (1-based)
nbdry_ecco = 300

# upper and lower bounds for temp and salinity
tempUp = 10
tempLow = -3

saltUp = 34.8
saltLow = 33.2

# define grid bounds of lake vostock
vostock=array([545,325,583,347])

# call function that generates the out file
roms_ini_ecco(grd_file, theta_file, salt_file, out_file, Tcline, theta_s, theta_b, hc, N, nbdry_ecco,(tempLow,tempUp),(saltLow,saltUp),vostock)
# set user parameter and call main routine

grd_file = os.path.join(os.environ['projdir'],'data','preprocessing','processed','waom10_grd.nc')
out_file = os.path.join(os.environ['projdir'],'data','preprocessing','processed','waom10_ini.nc')

print('making ini file for grid file: ' + grd_file)

# Path to ECCO2 files for temperature and salinity 1th January 1995
salt_file = os.path.join(os.environ['projdir'],'data','preprocessing','external','ecco2','SALT.nc','SALT.1440x720x50.20070101.nc')
theta_file = os.path.join(os.environ['projdir'],'data','preprocessing','external','ecco2','THETA.nc','THETA.1440x720x50.20070101.nc')

# Grid parameters; check grid_file and *.in file to make sure these are correct
Tcline = 50
theta_s = 4.0
theta_b = 0.9
hc = 50
N = 31

# Northernmost index of ECCO2 grid to read (1-based)
nbdry_ecco = 300

# upper and lower bounds for temp and salinity
tempUp = 10
tempLow = -3

saltUp = 34.8
saltLow = 33.2

# define grid bounds of lake vostock
vostock=array([545,325,583,347])

# call function that generates the out file
roms_ini_ecco(grd_file, theta_file, salt_file, out_file, Tcline, theta_s, theta_b, hc, N, nbdry_ecco,(tempLow,tempUp),(saltLow,saltUp),vostock)



from netCDF4 import Dataset
from numpy import *
import os
import sys
from scipy.interpolate import NearestNDInterpolator, RegularGridInterpolator
sose_path = os.path.join(os.environ['projdir'],'data','preprocessing','external','sose')
sys.path.append(sose_path)
from mds import *
import scipy.io as sio

#load roms
print('loading data: roms grid, sose salt and theta, sose grid')
grd_file = os.path.join(os.environ['projdir'],'data','preprocessing','processed','waom10_small_grd.nc')
out_file = os.path.join(os.environ['projdir'],'data','preprocessing','processed','waom10_small_nudge.nc')
sose_path = os.path.join(os.environ['projdir'],'data','preprocessing','external','sose')
salt_path = os.path.join(sose_path,'SALT_mnthlyBar')
theta_path = os.path.join(sose_path,'THETA_mnthlyBar')
grid_path = os.path.join(sose_path,'grid.mat')

id = Dataset(grd_file,'r')
zice = id.variables['zice'][:,:]
mask_rho = id.variables["mask_rho"][:,:]
lat_roms = id.variables['lat_rho'][:,:]
lon_roms = id.variables['lon_rho'][:,:]
id.close()

salt_raw = rdmds(salt_path,itrs=100,rec=np.arange(24,36),returnmeta=False,lev=[0],fill_value=NaN)
theta_raw = rdmds(theta_path,itrs=100,rec=np.arange(24,36),returnmeta=False,lev=[0],fill_value=NaN)
sose_grid = sio.loadmat(grid_path)

print('prepare sose data for interpolation')
#apply sose mask to sose data
sose_mask_raw = sose_grid["maskCtrlC"]
sose_mask = tile(swapaxes(sose_mask_raw[:,:,0],0,1),(12,1,1))

salt_ma = ma.masked_where(sose_mask==0.0,salt_raw)
theta_ma = ma.masked_where(sose_mask==0.0,theta_raw)

# load lon and lat sose and change lon to -180 to 180
lon_sose_raw = sose_grid["XC"][:,0]
lon_sose_raw[lon_sose_raw>180] -=360
lat_sose = sose_grid["YC"][0,:]

#reorder lon so it's strictly ascending
lon_order = argsort(lon_sose_raw)
lon_sose_raw = lon_sose_raw[lon_order]

# sose doesnt wrap around, so copy beginning and end
lon_sose_tmp = zeros(size(lon_sose_raw)+2)
lon_sose_tmp[0] = lon_sose_raw[-1]-360
lon_sose_tmp[1:-1] = lon_sose_raw
lon_sose_tmp[-1] = lon_sose_raw[0]+360
lon_sose = lon_sose_tmp.copy()

#reorder and copy sose_data according to lon manipulations
def reorder_sose(data):
    sss_sose_raw = ma.copy(data)
    sss_sose_raw = sss_sose_raw[:,:,lon_order]
    sss_sose_tmp = ma.zeros((size(sss_sose_raw,0),size(sss_sose_raw,1),size(sss_sose_raw,2)+2))
    sss_sose_tmp[:,:,0] = sss_sose_raw[:,:,-1]-360
    sss_sose_tmp[:,:,1:-1] = sss_sose_raw
    sss_sose_tmp[:,:,-1] = sss_sose_raw[:,:,0]+360
    sss_sose = sss_sose_tmp.copy()
    return sss_sose

salt = reorder_sose(salt_ma)
theta = reorder_sose(theta_ma)

#interpolate sose to roms grid
print('interpolate sose to roms grid and fill in mask')
def NDinterp(data):
            
    valid_mask = ~np.isnan(data)
    coords = np.array(np.nonzero(valid_mask)).T
    values = data[valid_mask]

    it = NearestNDInterpolator(coords,values)

    filled = it(list(np.ndindex(data.shape))).reshape(data.shape)

    return filled
    

def sose2roms(sose_data):

    sss_interp = ma.zeros((12,size(lat_roms,0),size(lat_roms,1)))

    for month,A in enumerate(sose_data):

            print("processing month: ",month)
    
            # fill in land mask with nearest neighbours
            print("fill in land mask")
            A[A.mask]=np.nan
            B = NDinterp(A)
            
            #interpolate to roms grid
            print("interpolate to roms grid")
            interp_func = RegularGridInterpolator((lat_sose,lon_sose),A,bounds_error=False, method="nearest",fill_value=NaN)
            C = interp_func((lat_roms,lon_roms))

            #fill in far south region
            print("fill in far south")
            D = NDinterp(C)
            
            #D[D<bounds[0]]=bounds[0]
            #D[D>bounds[1]]=bounds[1]

            sss_interp[month] = D
            
    return sss_interp

salt_it = sose2roms(salt)
theta_it = sose2roms(theta)

print('set up dQdSST array and time array')
#set up surface net heat flux sensitivity to SST with dQdSST = -40 in takeshi melt season (Nov till Feb)
dQdSST=np.ones(np.shape(salt_it))*-40
dQdSST[:,zice<0.0]=0.0
dQdSST[:,mask_rho==0]=0.0
dQdSST[2:-2,lat_roms<=-55]=0.0

time_start = 365/12*0.5
time_step = 365/12
time = np.arange(time_start,365,time_step)

# Set up output file
num_lon = size(lon_roms, 1)
num_lat = size(lon_roms, 0)

print('Writing ' + out_file)
out_id = Dataset(out_file, 'w')
out_id.createDimension('xi_rho', num_lon)
out_id.createDimension('eta_rho', num_lat)
out_id.createDimension('sss_time', len(time))
out_id.createVariable('sss_time', 'f8', ('sss_time'))
out_id.createDimension('sst_time', len(time))
out_id.createVariable('sst_time', 'f8', ('sst_time'))
out_id.variables['sss_time'].long_name = 'time since initialization'
out_id.variables['sss_time'].units = 'days'
out_id.variables['sss_time'].cycle_length = 365.
out_id.variables['sst_time'].long_name = 'time since initialization'
out_id.variables['sst_time'].units = 'days'
out_id.variables['sst_time'].cycle_length = 365.
out_id.createVariable('SSS', 'f8', ('sss_time', 'eta_rho', 'xi_rho'))
out_id.variables['SSS'].long_name = 'surface salinity'
out_id.variables['SSS'].units = 'PSU'
out_id.createVariable('SST', 'f8', ('sst_time', 'eta_rho', 'xi_rho'))
out_id.variables['SST'].long_name = 'surface temperature'
out_id.variables['SST'].units = 'degree Celsius'
out_id.createVariable('dQdSST', 'f8', ('sst_time', 'eta_rho', 'xi_rho'))
out_id.variables['dQdSST'].long_name = 'surface net heat flux sensitivity to SST'
out_id.variables['dQdSST'].units = 'watt meter-2 Celsius-1'

out_id.variables['sss_time'][:] = time  
out_id.variables['sst_time'][:] = time        
out_id.variables['SSS'][:] = salt_it
out_id.variables['SST'][:] = theta_it
out_id.variables['dQdSST'][:] = dQdSST

out_id.close()

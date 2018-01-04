from netCDF4 import Dataset
import numpy as np
from scipy.spatial import KDTree
import scipy.io as sio

def read_atg(atg_path,site_id):
    mat_content = sio.loadmat(atg_path)
    atg_data = mat_content['atg']
    tide_data = {}
    for key in ['name','lat','lon','amp','Gphase']:
        tide_data[key] =np.squeeze(atg_data[key][0,0][site_id-1])
        tide_data['constit']=np.squeeze(atg_data['constit'][0,0][:])
    tide_data['name'] = tide_data['name'].strip()
    return tide_data

def extract_zeta(file_path,grid_path,target_lat,target_lon):

    print('extract zeta ...')
    
    #read roms values and grid data
    data = Dataset(file_path,'r')
    zeta_rho = data.variables['zeta'][:,:,:]
    data.close()

    grid = Dataset(grid_path,'r')
    lon_rho = grid.variables['lon_rho'][:,:]
    lat_rho = grid.variables['lat_rho'][:,:]
    mask_rho = grid.variables['mask_rho'][:,:]
    grid.close()

    lat_s = lat_rho[mask_rho==1]
    lon_s = lon_rho[mask_rho==1]
    zeta_s = zeta_rho[:,mask_rho==1]
    #define target point
    lat_t = target_lat
    lon_t = target_lon

    #build KDTree to perform nearest neighbour look up 
    point_list = np.column_stack((lat_s,lon_s))
    tree = KDTree(point_list)
    dist, ind = tree.query((lat_t,lon_t))
    print("target lat, lon: ",lat_t,lon_t)
    print("nearest neighbour lat, lon: ",lat_s[ind],lon_s[ind])

    #surface elevation timeseries at target location
    #zeta_s_flat = zeta_s.reshape((size(zeta_s,0),size(lon_s)))
    zeta_t = zeta_s[:,ind]

    return (zeta_t,ind)

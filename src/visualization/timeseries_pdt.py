
from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from rotate_vector_roms import *
from scipy.spatial import KDTree

# Calculate and plot timeseries of the Drake Passage transport during a
# ROMS simulation.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS history/averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def plot_dpt(file_path):
    
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Longitude of Drake Passage zonal slice (convert to ROMS bounds 0-360)
    lon0 = -67.0
    # Latitude bounds on Drake Passage zonal slice
    lat_min = -68.0
    lat_max = -54.5

    time = []
    dpt = []

    print( 'Reading grid')
    id = Dataset(file_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon_roms = id.variables['lon_rho'][:,:]
    lat_roms = id.variables['lat_rho'][:,:]
    angle = id.variables['angle'][:,:]
    mask = id.variables['mask_rho'][:,:]
    id.close()

    print( 'Reading data')
    id = Dataset(file_path, 'r')
    # Read time values and convert from seconds to years
    new_time = id.variables['ocean_time'][:]/(60*60*24*365.25)
    num_time = size(new_time)
    # Calculate time-dependent water column thickness: h + zice + zeta
    zeta = id.variables['zeta'][:,:,:]
    wct_tmp = tile(h, (num_time,1,1)) + tile(zice, (num_time,1,1)) + zeta
    # Read barotropic velocities in x-y space
    ubar_xy = id.variables['ubar'][:,:,:]
    vbar_xy = id.variables['vbar'][:,:,:]
    id.close()
    #print(shape(ubar),shape(vbar))

    print( 'Rotating velocity vector')
    # Rotate one time index at a time
    ubar_tmp=empty(shape(zeta))
    vbar_tmp=empty(shape(zeta))
    for t in range(num_time):        
        ubar_tmp[t], vbar_tmp[t] = rotate_vector_roms(ubar_xy[t,:,:], vbar_xy[t,:,:], angle)

    print('interpolate everything to dp coordinates')

    lat_DP=arange(lat_min,lat_max,0.10)
    lon_DP=lon0

    lat_s = lat_roms[mask==1]
    lon_s = lon_roms[mask==1]
    wct_s = wct_tmp[:,mask==1]
    ubar_s = ubar_tmp[:,mask==1]
    vbar_s = vbar_tmp[:,mask==1]
    #define target point
    lat_t = lat_DP
    lon_t = tile(lon_DP,(size(lat_t)))

    #build KDTree to perform nearest neighbour look up 
    point_list = column_stack((lat_s,lon_s))
    tree = KDTree(point_list)
    dist, ind = tree.query(column_stack((lat_t,lon_t)))

    wct_DP = wct_s[:,ind]
    ubar_DP = ubar_s[:,ind]
    vbar_DP = vbar_s[:,ind]

    # Calculate dy
    # First calculate latitude on edges of each cell
    middle_lat = 0.5*(lat_DP[:-1] + lat_DP[1:])
    s_bdry = 2*lat_DP[0] - middle_lat[0]
    n_bdry = 2*lat_DP[-1] - middle_lat[-1]
    lat_edges = zeros(size(lat_DP)+1)
    lat_edges[0] = s_bdry
    lat_edges[1:-1] = middle_lat
    lat_edges[-1] = n_bdry
    # Now calculate difference in latitude across each cell
    dlat_DP = lat_edges[1:] - lat_edges[:-1]
    # Convert to Cartesian space for dy in metres
    dy_DP = r*dlat_DP*deg2rad
    print("mean dy:", mean(dy_DP))

    for t in range(num_time):
        # Integrate ubar*wct*dy and convert to Sv
        dpt.append(sum(ubar_DP[t,:]*wct_DP[t,:]*dy_DP)*1e-6)

    print( 'Plotting')
    plot(new_time,dpt)
    xlabel('Years')
    ylabel('Drake Passage Transport (Sv)')
    grid(True)
    show()

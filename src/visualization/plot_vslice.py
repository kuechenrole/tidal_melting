# Plot the terrain-following vertical levels through a line given in
# roms grid coordinates [i_min,j_min] [i_max,j_max]
# This is a good way to test out different choices
# for vertical stretching parameters.
# Input:
# grid_path = path to ROMS grid file
# lon0 = longitude to interpolate to (-180 to 180)
# depth_min = deepest depth to plot (negative, in metres)
# Vstretching, theta_s, theta_b, hc, N = vertical stretching parameters (see
# *.in configuration file if unsure)

from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *


def plot_vslice (file_path, variable, tstep, depth_min, depth_max, i_min, j_min, i_max, j_max, Vstretching, theta_s, theta_b, hc, N, tstop=None, vmin=None, vmax=None,title_str=None):

    #read grid and variable at timestep
    id = Dataset(file_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    mask = id.variables['mask_rho'][:,:]
    var = id.variables[variable]
    if tstop==None:
        data = var[tstep]
    else:
        data = mean(var[tstep:tstop],axis=0)
    if hasattr(var, 'units'):
        unit=var.units
    else:
        unit='nondimensional'
    name = var.long_name
    id.close()

    #calc 3D field of depth values
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, None, Vstretching)

    #get a 3D land mask
    mask_3d = tile(mask,(N,1,1))

    #mask out land
    data_3d = ma.masked_where(mask_3d==0,data)

    #simple index array for later plotting
    idx_4d=np.indices(shape(data_3d),int)

    #extract values along this line (from stackoverflow)
    #make a line with num points
    num = int(sqrt(square(i_max-i_min)+square(j_max-j_min)))+1
    x,y=np.linspace(i_min,i_max,num),np.linspace(j_min,j_max,num)

    #extract values of 3d arrays along this line 
    z_2d =z_3d[:,y.astype(np.int),x.astype(np.int)]
    idx_3d =idx_4d[:,:,y.astype(np.int),x.astype(np.int)]
    data_2d =data_3d[:,y.astype(np.int),x.astype(np.int)]

    #convert index value to distance
    dist_2d=sqrt(square(idx_3d[1])+square(idx_3d[2]))

    #contour levels
    #lev=range(1,N)

    #Plot
    fig=figure(figsize=(18,6))
    if (vmin!=None and vmax!=None):
        pcolormesh(dist_2d,z_2d,data_2d,vmin=vmin,vmax=vmax)
    else:
        pcolormesh(dist_2d,z_2d,data_2d)
    colorbar()
    title(title_str)
    #title(name +" ("+ unit +") along line \n"+str([i_min, j_min])+" to "+str([i_max, j_max])+" (grid coords [i,j])")
    xlabel('Distance (km)')
    ylabel('Depth (m)')
    ylim([depth_min,depth_max])
    show()

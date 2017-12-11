from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_roms import *
import cmocean

# Make a circumpolar Antarctic plot of speed overlaid with velocity vectors at
# the given depth (surface, bottom, or vertically averaged).
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ocean history/averages file
# tstep = timestep in file_path to plot (1-indexed)
# depth_key = integer flag indicating whether to plot the surface velocity (1),
#             the bottom velocity (2), or vertically averaged velocity (3)
# save = optional boolean flag indicating that the plot should be saved to a
#        file rather than displayed on the screen
# fig_name = if save=True, filename for figure
def plot_avg (file_path, tstart, tstop):

    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Side length of blocks to average vectors over (can't plot vector at every
    # single point or the plot will be way too crowded)
    block = 15

    # Read angle from grid file
    #grid_id = Dataset(grid_path, 'r')
    #angle = grid_id.variables['angle'][j_min:j_max,i_min:i_max]
    #grid_id.close()
    # Read grid and velocity data
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    zeta = mean(id.variables['zeta'][tstart:tstop+1,:,:],axis=0)

    # Vertically averaged u and v
    u = mean(id.variables['ubar'][tstart:tstop+1,:,:],axis=0)
    v = mean(id.variables['vbar'][tstart:tstop+1,:,:],axis=0)
    id.close()
    
    angle = zeros(shape(lon))
    u_rho,v_rho = rotate_vector_roms(u, v, angle)

    speed = sqrt(square(u_rho)+square(v_rho))

    print('initialize and fill up the arrays')
    numy = size(lon,0) #530
    numx = size(lon,1) #630
    x = arange(numx)
    y = arange(numy)
    xmesh,ymesh = meshgrid(x,y)
    #print(numx,numy,x,y)

    # Average x, y, u_circ, and v_circ over block x block intervals
    # Calculate number of blocks
    size0 = int(ceil(numy/float(block)))
    size1 = int(ceil(numx/float(block)))
    # Set up arrays for averaged fields
    x_block = ma.empty([size0, size1])
    y_block = ma.empty([size0, size1])
    u_block = ma.empty([size0, size1])
    v_block = ma.empty([size0, size1])
    # Set up arrays containing boundary indices
    posn0 = list(arange(0, numy, block))
    posn0.append(numy)
    posn1 = list(arange(0, numx, block))
    posn1.append(numx)
    #print(posn0,posn1)
    # Double loop to average each block (can't find a more efficient way to do
    # this)
    for j in arange(size0):
        for i in arange(size1):
            start0 = posn0[j]
            end0 = posn0[j+1]
            start1 = posn1[i]
            end1 = posn1[i+1]
            x_block[j,i] = mean(xmesh[start0:end0, start1:end1])
            y_block[j,i] = mean(ymesh[start0:end0, start1:end1])
            u_block[j,i] = mean(u_rho[start0:end0, start1:end1])
            v_block[j,i] = mean(v_rho[start0:end0, start1:end1])

    # Make the plot
    fig,(ax0,ax1) = subplots(2,figsize=(15,15))

    speedP = ax0.pcolormesh(xmesh,ymesh,speed*100, vmin=0,vmax=30, cmap=cmocean.cm.speed)
    colorbar(speedP,ax=ax0)
    # Add vectors for each block
    quiverP = ax0.quiver(x_block, y_block, u_block, v_block,pivot="mid", color='black',units="width")
    quiverkey(quiverP, 0.7, 0.99, 0.2, r'$20 \frac{cm}{s}$', labelpos='E',
                       coordinates='figure')
    ax0.set_title('Mean barotropic velocity (cm/s)', fontsize=16)
    ax0.set_aspect('equal')
    ax0.axis('off')

    sshP = ax1.pcolormesh(zeta,vmin=-2,vmax=2,cmap=cm.bwr)
    colorbar(sshP,ax=ax1)
    ax1.set_title("Mean sea surface height [m]", fontsize=16)
    ax1.set_aspect("equal")
    ax1.axis("off")

    tight_layout()


    show()

def plot_widget (file_path, tstart,tstop):
    
    from ipywidgets import interact
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Side length of blocks to average vectors over (can't plot vector at every
    # single point or the plot will be way too crowded)
    block = 15
    
    print('read in the data')
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    zeta = id.variables['zeta'][tstart:tstop+1,:,:]

    # Vertically averaged u and v
    u = id.variables['ubar'][tstart:tstop+1,:,:]
    v = id.variables['vbar'][tstart:tstop+1,:,:]
    id.close()
    
    print('initialize and fill up the arrays')
    numt = size(u,0)
    numy = size(lon,0) #530
    numx = size(lon,1) #630
    
    u_rho = ma.empty([numt,numy,numx])
    v_rho = ma.empty([numt,numy,numx])
    speed = ma.empty([numt,numy,numx])
    
    angle = zeros(shape(lon))
    
    x = arange(numx)
    y = arange(numy)
    xmesh,ymesh = meshgrid(x,y)
    #print(numx,numy,x,y)

    # Average x, y, u_circ, and v_circ over block x block intervals
    # Calculate number of blocks
    sizet = size(u,0)
    size0 = int(ceil(numy/float(block)))
    size1 = int(ceil(numx/float(block)))
    # Set up arrays for averaged fields
    x_block = ma.empty([size0, size1])
    y_block = ma.empty([size0, size1])
    u_block = ma.empty([sizet,size0, size1])
    v_block = ma.empty([sizet,size0, size1])
    # Set up arrays containing boundary tices
    posn0 = list(arange(0, numy, block))
    posn0.append(numy)

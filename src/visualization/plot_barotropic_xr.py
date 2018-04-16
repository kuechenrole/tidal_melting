import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cmo 
import os
import sys
import xarray as xr

# add the 'src/' directory as one where we can import modules
src_dir = os.path.join(os.environ.get('projdir'),'src')
sys.path.append(src_dir)

from features.rotate_vector_roms import rotate_vector_roms

def plot_uv(ds,grd,block=15):
    
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = np.pi/180
    # Side length of blocks to average vectors over (can't plot vector at every
    # single point or the plot will be way too crowded)

    lon = grd.lon_rho.values
    lat = grd.lat_rho.values

    u = ds.ubar.values
    v = ds.vbar.values

    angle = np.zeros(np.shape(lon))
    u_rho,v_rho = rotate_vector_roms(u, v, angle)

    speed = np.sqrt(np.square(u_rho)+np.square(v_rho))

    #print('initialize and fill up the arrays')
    numy = np.size(lon,0) #530
    numx = np.size(lon,1) #630
    x = np.arange(numx)
    y = np.arange(numy)
    xmesh,ymesh = np.meshgrid(x,y)
    #print(numx,numy,x,y)

    # Average x, y, u_circ, and v_circ over block x block intervals
    # Calculate number of blocks
    size0 = int(np.ceil(numy/float(block)))
    size1 = int(np.ceil(numx/float(block)))

    # Set up arrays for averaged fields
    x_block = np.ma.empty([size0, size1])
    y_block = np.ma.empty([size0, size1])
    u_block = np.ma.empty([size0, size1])
    v_block = np.ma.empty([size0, size1])
    # Set up arrays containing boundary indices
    posn0 = list(np.arange(0, numy, block))
    posn0.append(numy)
    posn1 = list(np.arange(0, numx, block))
    posn1.append(numx)
    #print(posn0,posn1)
    # Double loop to average each block (can't find a more efficient way to do
    # this)
    for j in np.arange(size0):
        for i in np.arange(size1):
            start0 = posn0[j]
            end0 = posn0[j+1]
            start1 = posn1[i]
            end1 = posn1[i+1]
            x_block[j,i] = np.mean(xmesh[start0:end0, start1:end1])
            y_block[j,i] = np.mean(ymesh[start0:end0, start1:end1])
            u_block[j,i] = np.mean(u_rho[start0:end0, start1:end1])
            v_block[j,i] = np.mean(v_rho[start0:end0, start1:end1])

    # Make the plot
    fig,ax0 = plt.subplots(1,figsize=(15,10))

    speedP = ax0.pcolormesh(xmesh,ymesh,speed*100, vmin=0,vmax=30, cmap=cmo.speed)
    plt.colorbar(speedP,ax=ax0)
    # Add vectors for each block
    quiverP = ax0.quiver(x_block, y_block, u_block, v_block,pivot="mid", color='black',units="width")
    plt.quiverkey(quiverP, 0.75, 0.99, 0.1, r'$10 \frac{cm}{s}$', labelpos='E',
                       coordinates='figure')
    ax0.set_title('Mean barotropic velocity (cm/s)', fontsize=16)
    ax0.set_aspect('equal')
    ax0.axis('off')

    plt.show()

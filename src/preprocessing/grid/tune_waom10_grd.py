import sys
import os
from shutil import copyfile
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
src_dir = os.path.join(os.environ.get('projdir'),'src')
sys.path.append(src_dir)
from features.uvp_masks import uvp_masks

grd_old = os.path.join(os.environ.get('intdir'),'waom10_MinDepth20m_rx10.3_grd.nc')
grd_new = os.path.join(os.environ.get('prodir'),'waom10_grd.nc')
copyfile(grd_old,grd_new)

# mask zice north
def mask_zice(grd_file,lat):
    id=Dataset(grd_new,'r')
    zice = id.variables['zice'][:,:]
    lat_rho = id.variables["lat_rho"][:,:]
    id.close()

    mask_60 = np.where(lat_rho>=lat,0.0,1.0)
    zice_new=zice*mask_60

    id = Dataset(grd_new,'a')
    id.variables["zice"][:,:]=zice_new
    id.close()

def deepen_bathy(grid_file,min_depth):
    
    id = Dataset(grid_file,'a')
    h_old = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lat = id.variables['lat_rho'][:,:]
    mask_old = id.variables['mask_rho'][:,:]
    nbModif = 0
    nbModif_mask = 0
  # calculate watercolumn thickness and mask with land mask
    wc = h_old + zice
    h_new=h_old.copy()
    mask_new = mask_old.copy()
    for iEta in range(np.size(h_old,0)):
        for iXi in range(np.size(h_old,1)):
            if (mask_old[iEta,iXi]==1 and lat[iEta,iXi] <= -60.0 and wc[iEta,iXi]<min_depth):
                h_new[iEta,iXi] = h_old[iEta,iXi] + (min_depth - wc[iEta,iXi])
                nbModif += 1
                
            elif (mask_old[iEta,iXi]==1 and lat[iEta,iXi] > -60.0 and wc[iEta,iXi]<min_depth):
                mask_new[iEta,iXi]==0
                nbModif_mask+=1

    print('     nbModif=', nbModif)
    print('     nbModif_mask=',nbModif_mask)

    umask,vmask,pmask=uvp_masks(mask_new)

    id.variables['h'][:,:]= h_new
    id.variables['mask_rho']= mask_new
    id.variables['mask_u'][:,:]= umask
    id.variables['mask_v'][:,:]= vmask
    id.variables['mask_psi'][:,:]= pmask
        
    id.close()

def mask_box(grid_file,box):
    
    id = Dataset(grid_file,'a')
    zice = id.variables['zice'][:,:]
    mask_old = id.variables['mask_rho'][:,:]
    mask_new = mask_old.copy()
    nbModif = 0

    [imin,jmin,imax,jmax]=box
    imax=imax+1
    jmax=jmax+1

    for iEta in range(jmin,jmax):
        for iXi in range(imin,imax):
            if (mask_old[iEta,iXi]==1 and zice[iEta,iXi]<0.0):
                mask_new[iEta,iXi] = 0
                nbModif += 1

    print('     nbModif=', nbModif)


    umask,vmask,pmask=uvp_masks(mask_new)

    id.variables['mask_rho'][:,:]= mask_new
    id.variables['mask_u'][:,:]= umask
    id.variables['mask_v'][:,:]= vmask
    id.variables['mask_psi'][:,:]= pmask
        
    id.close()

mask_zice(grd_new,-60.0)
deepen_bathy(grd_new,50)

# boxes for small grid [imin,jmin,imax,jmax]

box1=np.array([481,117,520,137])
box2=np.array([190,590,210,607])
box3=np.array([167,532,167,532])
box4=np.array([170,558,170,558])
box5=np.array([695,357,699,362])

box_list = [box1,box2,box3,box4,box5]

# then do the masking
for box in box_list:
    mask_box(grd_new,box)

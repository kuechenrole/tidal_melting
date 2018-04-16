
# coding: utf-8

# In[11]:


import numpy as np 
import os
import sys
import xarray as xr
import scipy.io as sio
import matplotlib.pyplot as plt
import datetime

from dotenv import load_dotenv, find_dotenv

# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()
load_dotenv(dotenv_path)
src_dir = os.environ.get('srcdir')
sys.path.append(src_dir)

# always reload modules marked with "%aimport"
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '1')


from features.resample_roms import resample
from features.bathy_smoothing import smoothing_PlusMinus_rx0
from features.cartesian_grid_2d import haversine


# In[12]:


#run = os.environ.get('run')
run ='waom10'
mr = 10 #km
smooth = False
deepen = False


# In[13]:


#establish the grid with grid point distances of mr/2 in km
#we need double resolution to cover all of the staggered grid points (we subset to rho, psi, u, v points later)
#we need an extra line of u and v points at first to calculate all dx and dy on rho points

x,y = np.meshgrid(np.arange(-4300,4300+mr/2,mr/2),np.arange(-3700,3600+mr/2,mr/2))


# In[14]:


#load south polar stereographic projection to convert from grid point distance in m to lat/lon and back
from mpl_toolkits.basemap import Basemap
m = Basemap(projection='spstere',lon_0=0,boundinglat=-50,lat_ts=-71)

#get lat/lon coordinates at all grid points by shifting the grid to the lower left corner of the map
lon,lat=m(x*1000+m.urcrnrx/2,y*1000+m.urcrnry/2,inverse=True)

#calculate curvilinear coordinate distances at rho points
dx = haversine(lon[1::2,0:-2:2],lat[1::2,0:-2:2],lon[1::2,2::2],lat[1::2,2::2])
dy = haversine(lon[0:-2:2,1::2],lat[0:-2:2,1::2],lon[2::2,1::2],lat[2::2,1::2])

#calculate curvilinear coordinate metrices
pm = 1.0/dx
pn = 1.0/dy
 
dndx = np.empty_like(pm)
dmde = np.empty_like(pn)

dndx[:,1:-1] = 0.5*(pn[:,2:] - pn[:,:-2])
dmde[1:-1,:] = 0.5*(pm[2:,:] - pm[:-2,:])

dndx[:,0]  = 2*dndx[:,1]  - dndx[:,2]
dndx[:,-1] = 2*dndx[:,-2] - dndx[:,-3]
dmde[0,:]  = 2*dmde[1,:]  - dmde[2,:]
dmde[-1,:] = 2*dmde[-2,:] - dmde[-3,:]

#subset lat and lon at rho, psi, u and v points
lon_rho = lon[1::2,1::2]
lat_rho = lat[1::2,1::2]

lon_psi = lon[2:-1:2,2:-1:2]
lat_psi = lat[2:-1:2,2:-1:2]

lon_u = lon[1::2,2:-1:2]
lat_u = lat[1::2,2:-1:2]

lon_v = lon[2:-1:2,1::2]
lat_v = lat[2:-1:2,1::2]


# In[15]:


#load rtopo bed and ice topography and resample to rho points
rtopo_path = os.path.join(os.environ.get('extdir'),'rtopo','RTopo-2.0.1_30sec_*_S30.nc')
rtopo = xr.open_mfdataset(rtopo_path,data_vars='minimal')#.sel(latdim=np.arange(0,7501,50),londim=np.arange(0,43201,100))

rt_lon,rt_lat = np.meshgrid(rtopo.lon.values,rtopo.lat.values)
bed_raw = resample(rt_lon,rt_lat,lon_rho,lat_rho,rtopo.bedrock_topography.values)
ice_raw = resample(rt_lon,rt_lat,lon_rho,lat_rho,rtopo.ice_base_topography.values)

#make a copy of the raw bathymetry
bed = bed_raw.copy()
ice = ice_raw.copy()


# In[16]:


#set bed minimum depth to 10 cm
bed[bed>-0.1]= -0.1
#set ice draft at these places to zero 
ice[bed>0.1] = 0.0
#set ice mountains to zero
ice[ice>0]= 0.0

#set water column thickness to a small positive value (ROMS don't like when bed = ice draft)
wct = (-(bed-ice)).copy()
ice[wct==0] = bed[wct==0] + 0.1


# In[17]:


#generate a land/ocean mask depending on water column thickness
#(distance between ice and bed or sea surface and bed)
#wct = (-(bed-ice)).copy()
mask = np.ones_like(wct)
mask[wct<20] = 0


# In[18]:


if smooth:
    #if smoothing is activated smooth wct and bed and set ice draft as bed + wct
    mask = np.ones_like(wct)
    mask[wct<=0.1] = 0
    dA = 1.0/(pn*pm) 
    bed, HmodifVal, ValueFct = smoothing_PlusMinus_rx0(mask,bed,0.3,dA)
    wct, HmodifVal, ValueFct = smoothing_PlusMinus_rx0(mask,wct,0.3,dA)
    ice = bed + wct
    
    #update the minimum wct points as before
    bed[bed>-0.1]= -0.1
    ice[bed>0.1] = 0.0
    ice[ice>0]= 0.0
    wct = (-(bed-ice)).copy()
    ice[wct==0] = bed[wct==0] + 0.1
    
    #update the mask
    wct = (-(bed-ice)).copy()
    mask = np.ones_like(wct)
    mask[wct<20] = 0
    
#if deepening is activated, deepen the bed to a minimum water column thickness of 50m 
if deepen:
    shallow = (wct<50)&(wct>=20)
    bed[shallow] = ice[shallow]-50.0


# In[19]:


#set spherical flag to 1, since we're creating a curvilinear spherical grid
spherical_da = xr.DataArray(int(1),name='spherical',attrs={'flag_meanings': 'Cartesian spherical',
 'flag_values': np.array([0, 1], dtype=int),
 'long_name': 'grid type logical switch'})

xl = mr*np.size(lat_rho,1)*1000
xl_da = xr.DataArray(xl,name='xl',attrs={'long_name': 'basin length in the XI-direction', 'units': 'meter'} )
el = mr*np.size(lon_rho,0)*1000
el_da = xr.DataArray(el,name='el',attrs={'long_name': 'basin length in the ETA-direction', 'units': 'meter'} )

angle = lon_rho/180.0*np.pi
angle_da = xr.DataArray(angle,name='angle',dims=['eta_rho','xi_rho'],attrs={'long_name': 'angle between XI-axis and EAST', 'units': 'radians'})

pn_da = xr.DataArray(pn,name="pn",dims=['eta_rho','xi_rho'],attrs={'long_name': 'curvilinear coordinate metric in ETA', 'units': 'meter-1'})
pm_da = xr.DataArray(pm,name='pm',dims=['eta_rho','xi_rho'],attrs={'long_name': 'curvilinear coordinate metric in XI', 'units': 'meter-1'})

dmde_da = xr.DataArray(dmde,name='dmde',dims=['eta_rho','xi_rho'],attrs={'long_name': 'ETA-derivative of inverse metric factor pm', 'units': 'meter'})
dndx_da = xr.DataArray(dndx,name='dndx',dims=['eta_rho','xi_rho'],attrs={'long_name': 'XI-derivative of inverse metric factor nm', 'units': 'meter'})

f = 2*7.29e-5*np.sin(lat_rho*np.pi/180)
f_da = xr.DataArray(f,name='f',dims=['eta_rho','xi_rho'],attrs={'long_name': 'Coriolis parameter at RHO-points', 'units': 'second-1'})

h_da = xr.DataArray(-bed,name='h',dims=['eta_rho','xi_rho'],attrs={'long_name': 'model bathymetry at RHO-points', 'units': 'meter'})
hraw_da = xr.DataArray(-bed_raw,name='hraw',dims=['eta_rho','xi_rho'],attrs={'long_name': 'Working bathymetry at RHO-points', 'units': 'meter'})

zice_da = xr.DataArray(ice,name='zice',dims=['eta_rho','xi_rho'],attrs={'long_name': 'model ice draft at RHO-points', 'units': 'meter'})

lon_rho_da = xr.DataArray(lon_rho,name='lon_rho',dims=['eta_rho','xi_rho'],attrs={'long_name': 'longitude of RHO-points',
 'standard_name': 'longitude',
 'units': 'degree_east'})
lat_rho_da = xr.DataArray(lat_rho,name='lat_rho',dims=['eta_rho','xi_rho'],attrs={'long_name': 'latitude of RHO-points',
 'standard_name': 'latitude',
 'units': 'degree_north'})
lon_psi_da = xr.DataArray(lon_psi,name='lon_psi',dims=['eta_psi','xi_psi'],attrs={'long_name': 'longitude of psi-points',
 'standard_name': 'longitude',
 'units': 'degree_east'})
lat_psi_da = xr.DataArray(lat_psi,name='lat_psi',dims=['eta_psi','xi_psi'],attrs={'long_name': 'latitude of psi-points',
 'standard_name': 'latitude',
 'units': 'degree_north'})
lon_u_da = xr.DataArray(lon_u,name='lon_u',dims=['eta_u','xi_u'],attrs={'long_name': 'longitude of u-points',
 'standard_name': 'longitude',
 'units': 'degree_east'})
lat_u_da = xr.DataArray(lat_u,name='lat_u',dims=['eta_u','xi_u'],attrs={'long_name': 'latitude of u-points',
 'standard_name': 'latitude',
 'units': 'degree_north'})
lon_v_da = xr.DataArray(lon_v,name='lon_v',dims=['eta_v','xi_v'],attrs={'long_name': 'longitude of v-points',
 'standard_name': 'longitude',
 'units': 'degree_east'})
lat_v_da = xr.DataArray(lat_v,name='lat_v',dims=['eta_v','xi_v'],attrs={'long_name': 'latitude of v-points',
 'standard_name': 'latitude',
 'units': 'degree_north'})

from features.mask_roms_uvp import uvp_masks

mask_rho = mask.copy()
mask_u,mask_v,mask_psi = uvp_masks(mask_rho)

mask_rho_da = xr.DataArray(mask_rho,name='mask_rho',dims=['eta_rho','xi_rho'],attrs={'flag_meanings': 'land water',
 'flag_values': np.array([ 0.,  1.]),
 'long_name': 'mask on RHO-points'})
mask_psi_da = xr.DataArray(mask_psi,name='mask_psi',dims=['eta_psi','xi_psi'],attrs={'flag_meanings': 'land water',
 'flag_values': np.array([ 0.,  1.]),
 'long_name': 'mask on psi-points'})
mask_u_da = xr.DataArray(mask_u,name='mask_u',dims=['eta_u','xi_u'],attrs={'flag_meanings': 'land water',
 'flag_values': np.array([ 0.,  1.]),
 'long_name': 'mask on u-points'})
mask_v_da = xr.DataArray(mask_v,name='mask_v',dims=['eta_v','xi_v'],attrs={'flag_meanings': 'land water',
 'flag_values': np.array([ 0.,  1.]),
 'long_name': 'mask on v-points'})

grd = xr.Dataset({'spherical':spherical_da,
                'xl':xl_da,
                'el':el_da,
                'angle':angle_da,
                'pm':pn_da,
                'pn':pn_da,
                'dndx':dndx_da,
                'dmde':dmde_da,
                'f':f_da,
                'h':h_da,
                'hraw':hraw_da,
                'zice':zice_da,
                'lon_rho':lon_rho_da,
                'lat_rho':lat_rho_da,
                'lon_psi':lon_psi_da,
                'lat_psi':lat_psi_da,
                'lon_u':lon_u_da,
                'lat_u':lat_u_da,
                'lon_v':lon_v_da,
                'lat_v':lat_v_da,
                'mask_rho':mask_rho_da,
                'mask_psi':mask_psi_da,
                'mask_u':mask_u_da,
                'mask_v':mask_v_da,},
               attrs={'history': 'GRID file using make_grid.py, smoothing='+str(smooth)+
                      ', deepening='+str(deepen)+', '+str(datetime.date.today()),
                      'type': 'ROMS grid file'})


# In[20]:


out_path = os.path.join(os.environ.get('prodir'),'waom'+str(mr)+'_grd.nc')
#out_path = '~/raijin/short/m68/oxr581/waom10_test/waom10_grd_smooth.nc'
grd.to_netcdf(out_path,unlimited_dims='bath')


# Below just left overs from development
#include masks psi u v in smoothing and deepeningx_earth = np.empty((rows,columns))
y_earth = np.empty((rows,columns))

x_earth[:,x_0ind] = 0
y_earth[y_0ind,:] = 0

for column in np.arange(columns-x_0ind-1):
    dx = haversine(lon[:,x_0ind+column],lat[:,x_0ind+column],lon[:,x_0ind+column+1],lat[:,x_0ind+column+1])
    x_earth[:,x_0ind+column+1] = x_earth[:,x_0ind+column]+dx
    
for column in np.arange(x_0ind):    
    dx = haversine(lon[:,x_0ind-column],lat[:,x_0ind-column],lon[:,x_0ind-column-1],lat[:,x_0ind-column-1])
    x_earth[:,x_0ind-column-1] = x_earth[:,x_0ind-column]-dist
    
for row in np.arange(rows-y_0ind-1):  
    dy = haversine(lon[y_0ind+row,:],lat[y_0ind+row,:],lon[y_0ind+row+1,:],lat[y_0ind+row+1,:])
    y_earth[y_0ind+row+1,:] = y_earth[y_0ind+row,:]+dy

for row in np.arange(y_0ind):
    dy = haversine(lon[y_0ind-row,:],lat[y_0ind-row,:],lon[y_0ind-row-1,:],lat[y_0ind-row-1,:])
    y_earth[y_0ind-row-1,:] = y_earth[y_0ind-row,:]-dy
    def calc_angle(lon_u,lat_u,lon_v,lat_v):
    
    DTOR = np.pi/180
    
    a1 = lat_u[1:,:] - lat_u[:-1,:]
    a2 = lon_u[1:,:] - lon_u[:-1,:]
    a2[a2<-180]+=360
    a2[a2>180]-=360
    
    a2 = a2*np.cos(0.5*DTOR*(lat_u[1:,:]+lat_u[:-1,:]))
    angle = np.arctan2(a1,a2)
    
    a1 = lat_v[:,:-1] - lat_v[:,1:]
    a2 = lon_v[:,:-1] - lon_v[:,1:]
    a2[a2<-180]+=360
    a2[a2>180]-=360
    
    a2 = a2*np.cos(0.5*DTOR*(lat_v[:,:-1] + lat_v[:,1:]))
    angle = 0.5*(angle + np.arctan2(a2,-a1))    
    
    return angle

angle = calc_angle(lon_u,lat_u,lon_v,lat_v)

print(np.shape(angle))
plt.close()
plt.contour(angle)
plt.colorbar()
plt.show()x_rho_da = xr.DataArray(x_rho*1000,name='x_rho',dims=['eta_rho','xi_rho'],attrs={'long_name': 'X-location of RHO-points', 'units': 'meter'})
y_rho_da = xr.DataArray(y_rho*1000,name='y_rho',dims=['eta_rho','xi_rho'],attrs={'long_name': 'Y-location of RHO-points', 'units': 'meter'})x_psi_da = xr.DataArray(x_psi*1000,name='x_psi',dims=['eta_rho','xi_rho'],attrs={'long_name': 'X-location of PSI-points', 'units': 'meter'})
y_psi_da = xr.DataArray(y_psi*1000,name='y_psi',dims=['eta_psi','xi_psi'],attrs={'long_name': 'Y-location of PSI-points', 'units': 'meter'})x_u_da = xr.DataArray(x_u*1000,name='x_u',dims=['eta_u','xi_u'],attrs={'long_name': 'X-location of U-points', 'units': 'meter'})
y_u_da = xr.DataArray(y_u*1000,name='y_u',dims=['eta_u','xi_u'],attrs={'long_name': 'Y-location of U-points', 'units': 'meter'})x_v_da = xr.DataArray(x_v*1000,name='x_v',dims=['eta_v','xi_v'],attrs={'long_name': 'X-location of V-points', 'units': 'meter'})
y_v_da = xr.DataArray(y_v*1000,name='y_v',dims=['eta_v','xi_v'],attrs={'long_name': 'Y-location of V-points', 'units': 'meter'})%matplotlib notebook
#set up mask_rho depending on water column thickness
wct = (bed+ice).copy()
mask_rho = np.ones_like(wct)
mask_rho[wct<20] = 0

def make_neibhour_mask(mask):
    neib_xi = np.empty_like(mask)
    neib_eta = np.empty_like(mask)

    neib_xi[:,1:-1] = mask[:,2:]+mask[:,:-2]
    neib_xi[:,0] = neib_xi[:,1]
    neib_xi[:,-1] = neib_xi[:,-2]

    neib_eta[1:-1,:] = mask[2:,:]+mask[:-2,:]
    neib_eta[0,:] = neib_eta[1,:]
    neib_eta[-1,:] = neib_eta[-2,:]

    return neib_eta + neib_xi#first get rid of any dubious ocean point nort of 60S
mask_sa = mask_rho.copy()
mask_sa[(lat_rho>-60)&(neib<=1)] = 0

#neib = make_neibhour_mask(mask_sa)
#mask_sa[(lat_rho>-60)&(neib<=3)] = 0

plt.close()
plt.pcolormesh(mask_sa)
plt.colorbar()
plt.show()intdirintdirplt.close()
plt.pcolormesh(neib)
plt.colorbar()
plt.show()#then, more carefully, find and mask ocean point around antarctica that are sorounded by land
Nmodif=1
while Nmodif > 0:
    mask_new = mask_tmp.copy()
    neib = make_neibhour_mask(mask_tmp)
    mask_new[(lat_rho<=-60)&(neib<=1)] = 0
    Nmodif = np.sum(mask_tmp)-np.sum(mask_new)
    print(Nmodif)
    mask_tmp = mask_new.copy()

mask_tmp[(lat_rho<=-60)&(neib<=2)] =0

Nmodif=1
while Nmodif > 0:
    mask_new = mask_tmp.copy()
    neib = make_neibhour_mask(mask_tmp)
    mask_new[(lat_rho<=-60)&(neib<=1)] = 0
    Nmodif = np.sum(mask_tmp)-np.sum(mask_new)
    print(Nmodif)
    mask_tmp = mask_new.copy()
plt.close()
plt.pcolormesh(mask_rho)
plt.colorbar()
plt.show()%matplotlib notebook
#include the amundsen fix from Milan 2017

#first load and resample to whole roms rho grid
import scipy.io as sio
amundsen_fix_path = os.path.join(os.environ.get('extdir'),'rtopo','Bathymetry_ASE_Millan_et_al_2017.nc')
amu_fix = xr.open_dataset(amundsen_fix_path)
amundsen_fix_latlon_path = os.path.join(os.environ.get('extdir'),'rtopo','lon_lat_romain_grid.mat')
amu_fix_latlon = sio.loadmat(amundsen_fix_latlon_path)

amu_lon = amu_fix_latlon['lone']
amu_lat = amu_fix_latlon['late']
amu_lon[amu_lon>180]-=360.0
amu = resample(amu_lon,amu_lat,lon_rho,lat_rho,amu_fix.BED.values)

#then lookup indices of coorner coordinates and subset the working bathymetry just there 
from scipy.spatial import KDTree

points = np.column_stack((lon_rho.flatten(),lat_rho.flatten()))

tree = KDTree(points)

target = np.column_stack((amu_lon[0,0],amu_lat[0,0]))
llcDist, llcInd = tree.query(target)
llcInd = np.squeeze(llcInd)

target = np.column_stack((amu_lon[-1,-1],amu_lat[-1,-1]))
urcDist, urcInd = tree.query(target)
urcInd = np.squeeze(urcInd)

#establish grid coordinates at rho points only (needed for the amundsen fix)
xi_rho = np.arange((np.size(x,1)-1)/2,dtype=int)
eta_rho = np.arange((np.size(x,0)-1)/2,dtype=int)

xi_2d,eta_2d = np.meshgrid(xi_rho,eta_rho)
xiEta = np.column_stack((xi_2d.flatten(),eta_2d.flatten()))

xi_amu_ll,xi_amu_ur = xiEta[llcInd][0],xiEta[urcInd][0]
eta_amu_ll,eta_amu_ur = xiEta[llcInd][1],xiEta[urcInd][1]

bed[eta_amu_ll:eta_amu_ur+1,xi_amu_ll:xi_amu_ur+1] = amu[eta_amu_ll:eta_amu_ur+1,xi_amu_ll:xi_amu_ur+1]
bed=bed*-1#have a look at the amundsen patch
plt.close()
plt.pcolormesh(bed[eta_amu_ll-20:eta_amu_ur+21,xi_amu_ll-20:xi_amu_ur+21])
plt.colorbar()
plt.show()
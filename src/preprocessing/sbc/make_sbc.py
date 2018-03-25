
# coding: utf-8

# In[7]:


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

get_ipython().run_line_magic('aimport', 'features.resample')
from features.resample import resample
from features.grid_ttide import NDinterp
from features.log_progress import log_progress

run = os.environ.get('run')
#run ='waom10'


# In[8]:


#read in tamura land mask
T_mask_path = os.path.join(os.environ.get('extdir'),'tamura','EASE_landmask_H.data')
with open(T_mask_path,'rb') as fid:
    T_mask = np.fromfile(fid,count=(721*721),dtype='float32').reshape((721,721))
    T_mask = np.flipud(T_mask)


# In[9]:


#get tamura lat lon coordinates
T_lat_lon_path = os.path.join(os.environ.get('extdir'),'tamura','latlon.data')
with open(T_lat_lon_path,'rb') as fid:
    T_lat_lon = np.fromfile(fid,count=(721*721*2),dtype='float32').reshape((2,721,721))
T_lat,T_lon = (T_lat_lon[0],T_lat_lon[1])
T_lat = np.flipud(T_lat)
T_lon = np.flipud(T_lon)


# In[10]:


#read in era interim winds and resample from twice daily to daily
era_path = os.path.join(os.environ.get('extdir'),'era_interim','ERA_Interim_1992_2011.2daily.*winds.nc')
era_ds = xr.open_mfdataset(era_path,data_vars='minimal').sel(time='2007',latitude=slice(-30,-90)).resample(time='D').mean()


# In[11]:


#get era coordinates
era_lat = era_ds.latitude.values
era_lon = era_ds.longitude.values
era_lon[era_lon>180]-=360.0
E_lon,E_lat = np.meshgrid(era_lon,era_lat)


# In[12]:


#get roms grid
R_grid_path = os.path.join(os.environ.get('prodir'),run+'_grd.nc')
R_grid = xr.open_dataset(R_grid_path)
R_lon = R_grid.lon_rho.values
R_lat = R_grid.lat_rho.values
R_angle = R_grid.angle.values
R_ulon = R_grid.lon_u.values
R_vlon = R_grid.lon_v.values
R_ulat = R_grid.lat_u.values
R_vlat = R_grid.lat_v.values


# In[13]:


month = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
daysPerMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
#month = ['jan','feb']
#daysPerMonth = [1,2]
dayOfYear = 1

for month,days in zip(month,daysPerMonth):
    
    print('Processing month: ',month,'with days: ',days)
    
    daysOfYear = np.arange(dayOfYear,dayOfYear+days,dtype=int)
    
    print('Containing days of year: ',daysOfYear)

    # preparing empty dataset
    ds = xr.Dataset({'shflux':(['shf_time','eta_rho','xi_rho'], np.empty((days,R_grid.eta_rho.size,R_grid.xi_rho.size))),
                     'swflux':(['swf_time','eta_rho','xi_rho'], np.empty((days,R_grid.eta_rho.size,R_grid.xi_rho.size))),
                     'sustr':(['sms_time','eta_u','xi_u'], np.empty((days,R_grid.eta_u.size,R_grid.xi_u.size))),
                     'svstr':(['sms_time','eta_v','xi_v'], np.empty((days,R_grid.eta_v.size,R_grid.xi_v.size)))},
                   coords={'shf_time':(['shf_time'],daysOfYear),
                           'swf_time':(['swf_time'],daysOfYear),
                           'sms_time':(['sms_time'],daysOfYear)})
   
    #open Tamura month flux data 
    T_data_path = os.path.join(os.environ.get('extdir'),'tamura','TSDM2hb_2007_'+month+'.data')
    with open(T_data_path,'rb') as fid:
        T_data = np.swapaxes(np.fromfile(fid,count = days*6*721*721 ,dtype='float32').reshape(days,6,721,721),0,1)
    
    #looping over the days with running day-of-the-year and day-of-the-month index
    for Eidx,Tidx in zip(log_progress(daysOfYear-1,name='days'),np.arange(days)):
        
        #read in Tamura heat and fresh water flux and turn in right position
        shflux_tmp = np.flipud(T_data[0,Tidx])
        ssflux_tmp = np.flipud(T_data[2,Tidx])
        
        #fill in tamuar mask for later resampling
        shflux_tmp[T_mask==0] = np.nan
        shflux_tmp = NDinterp(shflux_tmp)
        
        ssflux_tmp[T_mask==0] = np.nan
        ssflux_tmp = NDinterp(ssflux_tmp)

        #resample to roms grid points
        shflux_tmp = resample(T_lon,T_lat,R_lon,R_lat,shflux_tmp)
        ssflux_tmp = resample(T_lon,T_lat,R_lon,R_lat,ssflux_tmp)
        
        #correct large summer heat flux values and save to dataset
        shflux_tmp[shflux_tmp > 0.0]*=0.5
        
        ds.shflux[Tidx] = shflux_tmp
        
        #convert to freshwater flux with convention positive up 'swf (E-P)',
        #that means a positive freshwater flux value results in positive salt flux value
        #and save to dataset
        refSalt = 34.4

        ds.swflux[Tidx] = ssflux_tmp/refSalt*100
        
        #resample era-interim winds to roms grid
        uwnd = resample(E_lon,E_lat,R_lon,R_lat,era_ds.u10[Eidx].values)
        vwnd = resample(E_lon,E_lat,R_lon,R_lat,era_ds.v10[Eidx].values)
        
        #rotate wind directions to roms grid
        uv = (uwnd+1j*vwnd)*np.exp(1j*-R_angle)
        uwnd = uv.real
        vwnd = uv.imag
        
        #convert to stress
        signu = np.sign(uwnd)
        signv = np.sign(vwnd)

        rhoAir = 1.3
        Cd = 1.4e-3

        taux = rhoAir*Cd*np.square(uwnd)*signu
        tauy = rhoAir*Cd*np.square(vwnd)*signv
        
        #resample to roms u and v grid and save to dataset
        taux = resample(R_lon,R_lat,R_ulon,R_ulat,taux)
        tauy = resample(R_lon,R_lat,R_vlon,R_vlat,tauy)
        
        ds.sustr[Tidx]=taux
        ds.svstr[Tidx]=tauy
        
    #add attributes to data set and data arrays
    ds.attrs={'title':'waom surface heat/fresh water fluxes and wind stress',
                          'date':str(datetime.date.today()),
                          'tamura_file':T_data_path,
                          'era-interim file':era_path,
                          'grid file':R_grid_path,
                          'type':'ROMS forcing file'}
    ds.shflux.attrs = {'long_name': 'surface net heat flux', 'units': 'Watts meter-2'}
    ds.swflux.attrs = {'long_name': 'surface freshwater flux (E-P)',
                       'negative': 'net precipitation',
                       'positive': 'net evaporation',
                       'units': 'centimetre day-1'}
    ds.sustr.attrs = {'long_name': 'surface u-momentum stress', 'units': 'Newton meter-2'}
    ds.svstr.attrs = {'long_name': 'surface u-momentum stress', 'units': 'Newton meter-2'}
    ds.sms_time.attrs = {'cycle_length': days,'long_name': 'surface momentum stress time','units': 'day'}
    ds.shf_time.attrs = {'cycle_length': days, 'long_name': 'surface heat flux time', 'units': 'day'}
    ds.swf_time.attrs = {'cycle_length': days,'long_name': 'surface freshwater flux time','units': 'day'}
    
    #save month as netcdf file
    out_path = os.path.join(os.environ.get('intdir'),run+'_sbc_'+month+'.nc') 
    print("Saving month to "+out_path)
    ds.to_netcdf(out_path,'w')
    
    #update the day of the year value for next month
    dayOfYear += days


# In[19]:


import os
import xarray as xr
import glob

run_name = os.environ.get('run')


# concatenate all monthly data into a single file
monthly_data_paths = os.path.join(os.environ.get('intdir'),run+'_sbc_*.nc')
for file in glob.glob(monthly_data_paths):
    print('loading: '+file)

ds = xr.open_mfdataset(monthly_data_paths,concat_dim=None)
ds = ds.chunk(chunks={'shf_time':3,'swf_time':3,'sms_time':3})
for time in ['shf_time','swf_time','sms_time']:
    ds[time]['cycle_length'] = ds[time].cycle_length.size
year_data_path = os.path.join(os.environ.get('prodir'),run+'_sbc_2007.nc')

print('and saving to '+year_data_path)
ds.to_netcdf(year_data_path)


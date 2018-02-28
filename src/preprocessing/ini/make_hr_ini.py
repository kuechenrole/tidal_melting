
# coding: utf-8

# # Generate initial conditions file for high resolution experiment
# We want to initialize the high res run with the solution of the 10km rum. Therefore we have to map 2D and 3D velocities, SSH, salinity and temperature fields on the new grid.
# ## Preparation of ini file
# Load high and low resolution grid, low resolution history and low resolution ini file. 

# In[35]:


# get environment variables
import os
import sys
from dotenv import load_dotenv, find_dotenv

# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)

sys.path.append(os.environ.get('srcdir'))

# always reload modules marked with "%aimport"
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '1')

from features.resample import low_to_high
get_ipython().run_line_magic('aimport', 'features.resample')


# In[36]:


import xarray as xr
import matplotlib.pyplot as plt

hr_grd_path = os.path.join(os.environ.get("prodir"),'waom5_grd.nc')
hr_grd = xr.open_dataset(hr_grd_path)

lr_grd_path = os.path.join(os.environ.get("prodir"),'waom10_grd.nc')
lr_grd = xr.open_dataset(lr_grd_path)

lr_his_path =  os.path.join(os.environ.get("rawdir"),'waom10','ocean_his_ini_0007.nc')
lr_his = xr.open_dataset(lr_his_path).isel(ocean_time=0)

lr_ini_path = os.path.join(os.environ.get("prodir"),'waom10_ini.nc')
lr_ini = xr.open_dataset(lr_ini_path)


# Prepare high resolution ini file by dropping all horizontal grid dependent variables from the low resolution ini file. 

# In[37]:


hr_ini = lr_ini.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
hr_ini


# ## Interpolate low resolution variables on high resolution grid

# ### Function call for: zeta, ubar, vbar and u, v, temp, salt 
# Get 2D and 3D high resolution data and assign to prepared ini file.

# In[38]:


# interpolate 3D variables
hr_ini['u'] = low_to_high(lr_his.u,lr_grd,hr_grd,'u',3);
hr_ini['v'] = low_to_high(lr_his.v,lr_grd,hr_grd,'v',3);
hr_ini['temp'] = low_to_high(lr_his.temp,lr_grd,hr_grd,'rho',3)
hr_ini['salt'] = low_to_high(lr_his.salt,lr_grd,hr_grd,'rho',3)
# interpolate 2D variables
hr_ini['zeta'] = low_to_high(lr_his.zeta,lr_grd,hr_grd,'rho',2);
hr_ini['ubar'] = low_to_high(lr_his.ubar,lr_grd,hr_grd,'u',2); 
hr_ini['vbar'] = low_to_high(lr_his.vbar,lr_grd,hr_grd,'v',2);


# Assign time dimension to new data arrays with the ocean_time set to the one from the low resolution solution.

# In[21]:


import pandas as pd
hr_ini.coords['ocean_time'] = lr_his.ocean_time
hr_ini['ocean_time'] = pd.datetime(2007,1,1)
for var in ['zeta','ubar','vbar','u','v','temp','salt']:
    hr_ini[var] = hr_ini[var].expand_dims('ocean_time',0)


# ### Save new ini file
# Compare high resolution and low resolution ini file for consistency and save as netcdf file making sure that ocean_time is saved as unlimited dimension.

# In[22]:


hr_ini_path =  os.path.join(os.environ.get('prodir'),'waom5_ini.nc')
hr_ini.to_netcdf(hr_ini_path,unlimited_dims='ocean_time')


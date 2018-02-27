
# coding: utf-8

# # Generate initial conditions file for high resolution experiment
# We want to initialize the high res run with the solution of the 10km rum. Therefore we have to map 2D and 3D velocities, SSH, salinity and temperature fields on the new grid.
# ## Preparation of ini file
# Load high and low resolution grid, low resolution history and low resolution ini file. 

# In[3]:


# get environment variables
import os
from dotenv import load_dotenv, find_dotenv

# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)


# In[4]:


import xarray as xr
import numpy as np
import os
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

# In[5]:


hr_ini = lr_ini.drop(['u','v','ubar','vbar','salt','temp','zeta','ocean_time'])
hr_ini


# ## Interpolate low resolution variables on high resolution grid.
# ### Interpolation function
# Define function that takes low res data (as 2D or 3D data array), and high res grid and returns high resolution data array.

# In[12]:


from scipy import interpolate

def low_to_high(lr_da,lr_grd,hr_grd,gt,dim,fill_value=0.0):
    
    print('set up empty hr data array')
    if dim == 2:
    
        dummy = np.zeros(hr_grd['lon_'+gt].shape)
        x = hr_grd['xi_'+gt]
        y = hr_grd['eta_'+gt]
        hr_da = xr.DataArray(dummy,coords=[y,x],dims=['eta_'+gt,'xi_'+gt])
        
    elif dim == 3:
        
        N = lr_da.s_rho.size
        dummy = np.tile(np.zeros(hr_grd['lon_'+gt].shape),(N,1,1))
        x = hr_grd['xi_'+gt]
        y = hr_grd['eta_'+gt]
        z = lr_da['s_rho']
        hr_da = xr.DataArray(dummy,coords=[z,y,x],dims=['s_rho','eta_'+gt,'xi_'+gt])
    
    print('find index of shared first coordinate')
    # Find index of bottom left corner of low res data on high res grid
    ind = (hr_grd['x_'+gt].values == lr_grd['x_'+gt][0,0].values) & (hr_grd['y_'+gt].values == lr_grd['y_'+gt][0,0].values)
    eta0, xi0 = np.array(np.nonzero(ind)).squeeze()
    
    print('Fill in the mask of lr data')
    # Fill the mask of low resolution data with nearest neibghours and fill in known values on high res grid.
    if dim == 2:
        data = lr_da.values

        valid_mask = ~np.isnan(data)
        coords = np.array(np.nonzero(valid_mask)).T
        values = data[valid_mask]

        it = interpolate.NearestNDInterpolator(coords,values)

        filled = it(list(np.ndindex(data.shape))).reshape(data.shape)
        
        # Fill in known values on high res grid
        hr_da[eta0::2,xi0::2] = filled
        
    if dim == 3:
        
        for k in np.arange(N):
            
            print('processing depth level: ',k)
            data = lr_da[k].values

            valid_mask = ~np.isnan(data)
            coords = np.array(np.nonzero(valid_mask)).T
            values = data[valid_mask]

            it = interpolate.NearestNDInterpolator(coords,values)

            filled = it(list(np.ndindex(data.shape))).reshape(data.shape)
    
            # Fill in known values on high res grid
            hr_da[k,eta0::2,xi0::2] = filled
    
    # Now interpolate the intermediate grid points using x and y coordinates (in km from center point)
    print('interpolate at new intermediate cells on the hr grid')
    
    # Points we know the data
    x = hr_grd['x_'+gt][0,xi0::2].values
    y = hr_grd['y_'+gt][eta0::2,0].values
    
    # Define the target cells as meshgrid
    grid_x,grid_y = np.meshgrid(hr_grd['x_'+gt][0,:].values,hr_grd['y_'+gt][:,0].values)
    
    if dim == 2:
        # Data at these points
        values = hr_da[eta0::2,xi0::2].to_masked_array()

        # Define the interpolation function (fast regular grid interpolation)
        interp_func = interpolate.RegularGridInterpolator((y,x),values,bounds_error=False,fill_value=None)

        # Interpolate using linear interpolation
        interp = interp_func((grid_y[:,:],grid_x[:,:]))

        # Assign new data to data array
        hr_da[:,:] = interp
        
        # Fill with zeros where mask is present
        print('fill hr mask areas with fill value: ',fill_value)
        hr_da.values[hr_grd['mask_'+gt].values == 0] = fill_value
        
    if dim == 3:
        
        for k in np.arange(N):
            print('processing depth level: ',k)
            # Data at these points
            values = hr_da[k,eta0::2,xi0::2].to_masked_array()

            # Define the interpolation function (fast regular grid interpolation)
            interp_func = interpolate.RegularGridInterpolator((y,x),values,bounds_error=False,fill_value=None) 

            # Interpolate using linear interpolation
            interp = interp_func((grid_y[:,:],grid_x[:,:]))

            # Assign new data to data array
            hr_da[k,:,:] = interp

            # Fill with zeros where mask is present
            print('fill hr mask areas with fill value: ',fill_value)
            hr_da[k].values[hr_grd['mask_'+gt].values == 0] = fill_value
    
    return hr_da


# ### Function call for: zeta, ubar, vbar and u, v, temp, salt 
# Get 2D and 3D high resolution data and assign to prepared ini file.

# In[13]:


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


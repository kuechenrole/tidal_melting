# Generate initial conditions file for high resolution experiment
# We want to initialize the high res run with the solution of the 10km rum. Therefore we have to map 2D and 3D velocities, SSH, salinity and temperature fields on the new grid.
# ## Preparation of ini file
# Load high and low resolution grid, low resolution history and low resolution ini file. 

# In[35]:


# get environment variables
import os
from dotenv import load_dotenv, find_dotenv

# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)

hr_grd_path = os.path.join(os.environ.get("prodir"),'waom5_grd.nc')
lr_grd_path = os.path.join(os.environ.get("prodir"),'waom10_grd.nc')
lr_his_path =  os.path.join(os.environ.get("rawdir"),'waom10','ocean_his_ini_0007.nc')
lr_ini_path = os.path.join(os.environ.get("prodir"),'waom10_ini.nc')
hr_ini_path =  os.path.join(os.environ.get('prodir'),'waom5_ini.nc')

import make_hr_ini


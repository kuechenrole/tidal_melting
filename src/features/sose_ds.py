import os
import sys
import xarray as xr
import pandas as pd
import scipy.io as sio
import numpy as np
import features.mds as mds


def make_TS_ds(sose_dir=os.path.join(os.pardir,'data','external','sose'),records=None):
    '''Reads in SothernOceanStateEstimate Temperatures and Salinities and returns them in a Xarray dataset.'''
    
    #load grid data
    print("load grid")
    grid_path = os.path.join(sose_dir,'grid.mat')
    grid_raw = sio.loadmat(grid_path)
    
    XC = grid_raw["XC"]
    YC = grid_raw['YC']
    RC = grid_raw['RC'].squeeze()
    DRC = grid_raw['DRC'].squeeze()
    maskC = grid_raw['maskCtrlC']
    DXC = grid_raw['DXC']
    DYC = grid_raw['DYC']
    Depth = grid_raw['Depth']
    
    #load temperature data
    print("load temperature")
    temp_path = os.path.join(sose_dir,'THETA_mnthlyBar')
    temp_raw = mds.rdmds(temp_path,100,rec=records,fill_value=np.NaN)
    
    #load salt data
    print("load salt")
    salt_path = os.path.join(sose_dir,'SALT_mnthlyBar')
    salt_raw = mds.rdmds(salt_path,100,rec=records,fill_value=np.NaN)
    
    #define array of datetime range
    time_range = pd.period_range('2005-01',periods=len(temp_raw),freq='M')
    time_stamp = pd.Timestamp('2005-01')
    
    print("construct Xarray dataset")
    ds = xr.Dataset({'temperature':(['time','RC','YC','XC'],temp_raw),
                     'salinity':(['time','RC','YC','XC'],salt_raw),
                     'maskC':(['XC','YC','RC'],maskC),
                     'DRC':(['RC'],DRC),
                     'DXC':(['XC','YC'],DXC),
                     'DYC':(['XC','YC'],DYC),
                     'Depth':(['XC','YC'],Depth)},
                    coords={'latitude_center':(('XC','YC'),XC),
                           'longitude_center':(('XC','YC'),YC),
                           'depth_center':(('RC'),RC),
                           'time_range':(('time'),time_range),
                           'reference_time':time_stamp})
    print("done!")
    
    return ds
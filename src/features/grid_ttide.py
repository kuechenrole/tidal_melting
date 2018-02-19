import ttide as tt
import pandas as pd
from scipy.interpolate import NearestNDInterpolator
import xarray as xr
import numpy as np
from .log_progress import log_progress

def NDinterp(data):

    valid_mask = ~np.isnan(data)
    coords = np.array(np.nonzero(valid_mask)).T
    values = data[valid_mask]

    it = NearestNDInterpolator(coords,values)

    filled = it(list(np.ndindex(data.shape))).reshape(data.shape)

    return filled    


def grid_ttide(ds,grid_ds,const_list=['O1','M2'],res=50):
    
    ana_list = ['amp','amp_err','phase','phase_err']
    
    print('setting up the new fields ',ana_list,' for ',const_list)
    dummy = np.empty((ds.eta_rho.size,ds.xi_rho.size))
    dummy[:,:] = np.nan
    
    for const in const_list:
        for ana in ana_list:
            #print(const+'_'+ana)
            ds[const+'_'+ana]=(('eta_rho','xi_rho'),dummy.copy())
    
    
    stime = pd.to_datetime(ds.ocean_time[0].values)
     
    print("applying t_tide to every ",res,"th cell")
    xi_values = np.linspace(0,ds.xi_rho.size-1,res,dtype=int,endpoint=True)
    eta_values = np.linspace(0,ds.eta_rho.size-1,res,dtype=int,endpoint=True)
    
    for xi in log_progress(xi_values,name='xi'):
        
        for eta in eta_values:
            ds_sl = ds.isel(eta_rho=eta,xi_rho=xi)

            if ds_sl.zeta.isnull().values.any():
                for const in const_list:
                    for ana in ana_list:
                        ds[const+'_'+ana][eta,xi]=np.NaN
                
                
            else:
                signal = ds_sl.zeta.values
                latitude = ds_sl.lat_rho.values
                try:
                    ttide_out = tt.t_tide(signal,stime=stime,lat=latitude,out_style=None)

                    tt_ind_O1 = list(ttide_out['nameu']).index(b'O1  ')
                    tt_ind_M2 = list(ttide_out['nameu']).index(b'M2  ')
                    
                    tt_ind_list = [tt_ind_O1,tt_ind_M2]
                    
                    for const,tt_ind in zip(const_list,tt_ind_list):
                        for ana,tt_ana in zip(ana_list,ttide_out['tidecon'][tt_ind]):
                            ds[const+'_'+ana][eta,xi] = tt_ana

                except TypeError:
                    for const in const_list:
                        for ana in ana_list:
                            ds[const+'_'+ana][eta,xi]=np.NaN
                    
    print('interpolating intermediate cells and mask land')
    for con in const_list:
        for ana in ana_list:
            ds[con+'_'+ana].values = NDinterp(ds[con+'_'+ana].values)
            ds[con+'_'+ana] = ds[con+'_'+ana].where(grid_ds.mask_rho,0.0) 
      
        
    return ds

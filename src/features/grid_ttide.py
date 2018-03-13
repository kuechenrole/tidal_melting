import ttide as tt
import datetime
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


def grid_ttide(da,grid_ds,stime=datetime.datetime(2012,1,1),constit_list=['O1','M2'],res=50):
    
    ana_list = ['amp','amp_err','phase','phase_err']
    
    print('setting up the new fields ',ana_list,' for ',constit_list)
    dummy = np.empty((da.eta_rho.size,da.xi_rho.size))
    dummy[:,:] = np.nan
    
    for const in constit_list:
        for ana in ana_list:
            #print(const+'_'+ana)
            grid_ds[const+'_'+ana]=(('eta_rho','xi_rho'),dummy.copy())
     
    print("applying t_tide to every ",res,"th cell ..." )
    xi_values = np.linspace(da.xi_rho[0].values,da.xi_rho.size-1,res,dtype=int,endpoint=True)
    eta_values = np.linspace(da.eta_rho[0].values,da.eta_rho.size-1,res,dtype=int,endpoint=True)
    
    for xi in log_progress(xi_values,name='xi'):
        
        for eta in eta_values:
            da_sl = da.isel(eta_rho=eta,xi_rho=xi)

            if da_sl.isnull().values.any():
                for const in constit_list:
                    for ana in ana_list:
                        grid_ds[const+'_'+ana][eta,xi]=np.NaN
                
                
            else:
                signal = da_sl.values
                latitude = da_sl.lat_rho.values
                try:
                    ttide_out = tt.t_tide(signal,stime=stime,lat=latitude,out_style=None)
                    
                    tt_ind = {}
                    for const in constit_list:
                        tt_ind[const] = list(ttide_out['nameu']).index(str.encode(const+'  '))
                        
                        for ana,tt_ana in zip(ana_list,ttide_out['tidecon'][tt_ind[const]]):
                            grid_ds[const+'_'+ana][eta,xi] = tt_ana

                except TypeError:
                    for const in constit_list:
                        for ana in ana_list:
                            grid_ds[const+'_'+ana][eta,xi]=np.NaN
                    
    print('interpolating intermediate cells and mask land')
    for con in constit_list:
        for ana in ana_list:
            grid_ds[con+'_'+ana].values = NDinterp(grid_ds[con+'_'+ana].values)
            grid_ds[con+'_'+ana] = grid_ds[con+'_'+ana].where(grid_ds.mask_rho,0.0) 
      
        
    return grid_ds

import matplotlib.pyplot as plt

def plot_M2O1_diff(case_ds,case_str,ref_ds,ref_str,vmin=-0.10,vmax=0.10):
    
    plt.close('all')
    fig,axes = plt.subplots(ncols=2,nrows=2,figsize=(15,10))
    ax1,ax2,ax3,ax4 = axes.flatten()
    
    fig.suptitle('M2 and O1 height amplitude difference\n'+case_str+' - '+ref_str,fontsize=16)
     
    M2_diff = case_ds.M2_amp-ref_ds.tide_Eamp[0]
    O1_diff = case_ds.O1_amp-ref_ds.tide_Eamp[5]
    
    M2_diff_rel = (abs(case_ds.M2_amp-ref_ds.tide_Eamp[0])/ref_ds.tide_Eamp[0])*100
    O1_diff_rel = (abs(case_ds.O1_amp-ref_ds.tide_Eamp[5])/ref_ds.tide_Eamp[5])*100
    
    M2_diff.plot(ax=ax1,cmap=plt.cm.bwr,vmin=vmin,vmax=vmax)
    ax1.set_title('M2 ampl diff [m]')
    
    O1_diff.plot(ax=ax2,cmap=plt.cm.bwr,vmin=vmin,vmax=vmax)
    ax2.set_title('O1 ampl diff [m]')
    
    M2_diff_rel.fillna(0).plot(ax=ax3,vmin=0,vmax=100)
    ax3.set_title('M2 ampl relative diff [%]')
    
    O1_diff_rel.fillna(0).plot(ax=ax4,vmin=0,vmax=100)
    ax4.set_title('O1 ampl relative diff [%]')
    
    
    for ax in axes.flatten():
        ax.axis("off")
        ax.set_aspect('equal')
    
    plt.show()

def plot_M2O1_phase(case_ds,case_str,ref_ds,ref_str):
    plt.close('all')
    fig,axes = plt.subplots(ncols=2,nrows=2,figsize=(17,8))
    ax1,ax2,ax3,ax4 = axes.flatten()
    fig.suptitle('Comparison of M2 and O1 height phase\n '+case_str+' vs. '+ref_str,fontsize=16)

    case_ds.M2_phase.fillna(0).plot(ax=ax1)
    ax1.set_title(case_str +' [deg]')
    ax1.axis('off')

    ref_ds.tide_Ephase[0].plot(ax=ax2)
    ax2.set_title('TPXO M2 phase in deg')
    ax2.axis('off')
    
    case_ds.O1_phase.fillna(0).plot(ax=ax3)
    ax3.set_title(case_str+' [deg]')
    ax3.axis('off')

    ref_ds.tide_Ephase[5].plot(ax=ax4)
    ax4.set_title('TPXO O1 phase in deg')
    ax4.axis('off')

    for ax in axes.flatten():
        ax.set_aspect('equal')
        ax.axis("off")
        
    plt.show()

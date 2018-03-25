import xarray as xr

from .grid_ttide import grid_ttide,plot_amp,plot_phase
from .compare_atg import compare_atg,print_rmse

def compare_constit(case_amp,case_phase,case_str,ref_amp,ref_phase,ref_str,atg_rmse,wct,comp,constit):
    print_rmse(atg_rmse,[case_amp.name[:2]])
    plot_amp(case_amp,case_str,ref_amp,ref_str,comp,constit,wct)
    plot_phase(case_phase,case_str,ref_phase,ref_str,comp,constit)
    
def compare_tide(zeta,rds,stime,constits,stations,res,tpxods,case_str):

    print('tidal analysis at atg stations ...')
    comp,rmse = compare_atg(zeta,rds.mask_rho,stime=stime,constit_list=constits,station_list=stations,print_flag=False)
    
    print('tidal analysis on whole grid ...')
    rds = grid_ttide(zeta,rds,res=res)
    wct = rds.h+rds.zice

    print('write out ATG and TPXO comparison for constituent: ')
    for constit in constits:
        
        ind = ['M2','S2','N2','K2','K1','O1','P1','Q1'].index(constit)
        
        compare_constit(rds[constit+'_amp'],rds[constit+'_phase'],case_str,
                        tpxods.tide_Eamp[ind],tpxods.tide_Ephase[ind],'tpxo',
                       rmse,wct,comp,constit)
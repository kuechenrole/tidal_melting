from scipy.spatial import KDTree
import ttide as tt
import datetime
import numpy as np
import scipy.io as sio

def read_atg(atg_data,site_id,constit_list):

    site_data = {}
    for key in ['site_id','name','lat','lon','amp','Gphase','reclen','delta_t','meas_type','ref']:
        site_data[key] =np.squeeze(atg_data[key][0,0][site_id-1])
    site_data['constit']=np.squeeze(atg_data['constit'][0,0][:])
    site_data['name'] = site_data['name'].strip()

    cm2m = 1.0/100.0
    for const in constit_list:

        atg_con_ind = list(site_data['constit']).index(const)
        site_data[const]=np.array([site_data['amp'][atg_con_ind]*cm2m, site_data['Gphase'][atg_con_ind]])
        
    return site_data

def station_ttide(zeta_da,lat_t,lon_t,stime,constit_list):
    
    zeta_flat = zeta_da.stack(etaxi = ('eta_rho','xi_rho'))
    points = np.column_stack((zeta_flat.lat_rho.values,zeta_flat.lon_rho.values))
    tree = KDTree(points)
    
    target = np.column_stack((lat_t,lon_t))
    dist, ind = tree.query(target)
    
    tmp={}
    tmp['roms_signal'] = zeta_flat[:,ind].squeeze().values
    tmp['roms_ind'],tmp['dist_to ATG'] = ind,dist
    
    tmp['t_tide']=tt.t_tide(tmp['roms_signal'],dt=1,stime=stime,lat=zeta_flat.lat_rho[ind].values[0],out_style=None)


    for const in constit_list:
        tide_con_ind = list(tmp['t_tide']['nameu']).index(str.encode(const+'  ')) 
        tmp[const]=tmp['t_tide']['tidecon'][tide_con_ind]
        
    return tmp

def print_comparison(ttide_dict,atg_dict,constit_list):
      
        print(' Station: ',atg_dict['name'])
        print("Amp(amp_err)[m]:  atg      roms  || phase(phase_err)[deg]:  atg       roms")
        for con in constit_list:
            print(con,":             %0.2f"%atg_dict[con][0], "  %0.2f(%0.2f)"%(ttide_dict[con][0],ttide_dict[con][1]),\
                 "                      %0.2f"%atg_dict[con][1], "  %0.2f(%0.2f)"%(ttide_dict[con][2],ttide_dict[con][3])) 
            
def compare_tide(atg_mat_path,roms_zeta_da,stime=datetime.datetime(2007,1,1),constit_list = ['M2','O1'],station_list= [24,30]):

    mat_content = sio.loadmat(atg_mat_path)
    atg_data = mat_content['atg']
    
    for station in station_list:

        atg_dict = read_atg(atg_data,station,constit_list)
        lat = atg_dict['lat']
        lon = atg_dict['lon']
        ttide_dict = station_ttide(roms_zeta_da,lat,lon,stime,constit_list)

        print_comparison(ttide_dict,atg_dict,constit_list)
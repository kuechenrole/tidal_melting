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

def station_ttide(zeta_da,roms_mask_da,lat_t,lon_t,stime,constit_list):
    
    zeta_flat = zeta_da.stack(etaxi = ('eta_rho','xi_rho'))
    mask_flat = roms_mask_da.stack(etaxi = ('eta_rho','xi_rho'))
    
    lat_s = zeta_flat.lat_rho.values[mask_flat.values==True]
    lon_s = zeta_flat.lon_rho.values[mask_flat.values==True]
    zeta_s = zeta_flat.values[:,mask_flat.values==True]
    
    points = np.column_stack((lat_s,lon_s))
    tree = KDTree(points)
    
    target = np.column_stack((lat_t,lon_t))
    dist, ind = tree.query(target)
    
    tmp={}
    tmp['roms_signal'] = zeta_s[:,ind].squeeze()
    tmp['roms_ind'],tmp['dist_to ATG'] = ind,dist
    lat_r = lat_s[ind]
    lon_r = lon_s[ind]
    
    print('atg lat(lon): %.2f,%.2f'%(lat_t,lon_t))
    print('roms lat(lon): %.2f,%.2f'%(lat_r,lon_r))
    tmp['t_tide']=tt.t_tide(tmp['roms_signal'],dt=1,stime=stime,lat=lat_r,out_style=None)


    for const in constit_list:
        tide_con_ind = list(tmp['t_tide']['nameu']).index(str.encode(const+'  ')) 
        tmp[const]=tmp['t_tide']['tidecon'][tide_con_ind]
        
    return tmp

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def complex_rmse(predictions, targets):
    return np.sqrt(0.5*(((predictions - targets)*np.conjugate(predictions - targets)).real.mean()))

def calc_rmse(station_dict,constit_list):
    
    const_rmse={}

    for constit in constit_list:
        
        tt_amp_all = []
        atg_amp_all = []

        tt_phi_all =[]
        atg_phi_all =[]

        tt_z_all = []
        atg_z_all = []

        for station,data in station_dict.items():
            
            
            tt_amp = data['tt'][constit][0]
            atg_amp = data['atg'][constit][0]

            tt_phi = data['tt'][constit][2]
            atg_phi = data['atg'][constit][1]

            tt_amp_all.append(tt_amp)
            atg_amp_all.append(atg_amp)

            tt_phi_all.append(tt_phi)
            atg_phi_all.append(atg_phi)

            tt_z_all.append(tt_amp * np.exp(1j*tt_phi))
            atg_z_all.append(atg_amp * np.exp(1j*atg_phi))

        const_rmse[constit] = {}
        
        const_rmse[constit]['amp']=rmse(np.asarray(tt_amp_all),np.asarray(atg_amp_all))
        const_rmse[constit]['phase']=rmse(np.asarray(tt_phi_all),np.asarray(atg_phi_all))
        const_rmse[constit]['complex_amp']=complex_rmse(np.asarray(atg_z_all),np.asarray(tt_z_all))
        
    return const_rmse 

def print_station_dict(station_dict,constit_list):
    print("Amp(amp_err)[m]:  atg      roms    ||   phase(phase_err)[deg]:  atg       roms")
    for constit in constit_list:
        print(constit)
        for station_id,data in station_dict.items():  
            print('Station ',station_id,\
                  "      %0.2f"%data['atg'][constit][0], "  %0.2f(%0.2f)"%(data['tt'][constit][0],data['tt'][constit][1]),\
                  "                          %0.2f"%data['atg'][constit][1], "  %0.2f(%0.2f)"%(data['tt'][constit][2],data['tt'][constit][3])) 

def print_rmse(rmse_dict,constit_list):
    print('RMSE:  amp [m]    phase [deg]   complex_amp [m]')
    for constit,data in rmse_dict.items():
        print(constit,'       %.2f        %.2f         %.2f '%(data['amp'],data['phase'],data['complex_amp']))

def compare_tide(atg_mat_path,roms_zeta_da,roms_mask_da,stime=datetime.datetime(2007,1,1),constit_list = ['M2','O1'],station_list= [24,30]):

    print('stime = ',stime,' constits = ',constit_list,'stations = ',station_list)
    mat_content = sio.loadmat(atg_mat_path)
    atg_data = mat_content['atg']
    
    station_dict = {}
    
    for station in station_list:
        
        print('processing station ',station)
        station_dict[station] = {}

        atg_dict = read_atg(atg_data,station,constit_list)
        lat = atg_dict['lat']
        lon = atg_dict['lon']
        tt_dict = station_ttide(roms_zeta_da,roms_mask_da,lat,lon,stime,constit_list)

        #print_comparison(tt_dict,atg_dict,constit_list)
        
        station_dict[station]['atg'] = atg_dict
        station_dict[station]['tt'] = tt_dict
        
    print_station_dict(station_dict,constit_list)
    
   
        
    rmse_dict = calc_rmse(station_dict,constit_list)
    
    print_rmse(rmse_dict,constit_list)
    
    
        
    return station_dict,rmse_dict
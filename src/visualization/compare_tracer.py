#define function that generates a sose vs roms longitude transsect as monthly widget or annual meandefine function that generates a sose vs roms longitude transsect as monthly widget or annual mean
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import xarray as xr
import cmocean.cm as ocm
from ipywidgets import interact

def plot_lont(sds,rds,lon,lats,max_depth=None,tmin=-3,tmax=1,smin=33.8,smax=34.8,mean=True,month=None):
    
    if lon < 0.0:
        slon = lon + 360.0
        rlon = lon
    elif lon > 180:
        slon = lon
        rlon = lon - 360.0
    else:
        slon = lon
        rlon = lon
        
    print('lookup nearest neighbors from roms')
    rtemp_flat = rds.temp.stack(xieta=('xi_rho','eta_rho'))
    rsalt_flat = rds.salt.stack(xieta=('xi_rho','eta_rho')) 
    points = np.column_stack((rtemp_flat.lat_rho.values,rtemp_flat.lon_rho.values))
    tree = KDTree(points)
    
    print('lookup nearest neighbors from sose')
    lats_t = sds.sel(longitude=slon,latitude=lats,method='nearest').latitude.values
    lons_t = np.full(np.shape(lats_t),rlon,dtype=int)
    target = np.column_stack((lats_t,lons_t))
    dist, ind = tree.query(target)
    
    print('define axes and values from roms')
    rtemp_tr = rtemp_flat[:,:,ind]
    rsalt_tr = rsalt_flat[:,:,ind]
    x=rtemp_tr[0].lat_rho.fillna(0).values
    y=rtemp_tr[0].depth.fillna(0).values
    rtemp_val=rtemp_tr.to_masked_array()
    rsalt_val=rsalt_tr.to_masked_array()
    
    print('generate widget including sose plotting')
    plt.close()
    
    if mean==False:
        def plot(month):
            fig,axes = plt.subplots(2,2,figsize=(15,12))
            (ax1,ax2,ax3,ax4) = axes.flatten()

            sds.temperature.sel(longitude=slon,latitude=lats,method='nearest')[month].plot(ax=ax1,x='latitude',y='depth',vmin=tmin,vmax=tmax,cmap=ocm.thermal)
            ax1.set_title('SOSE')
            sds.salinity.sel(longitude=slon,latitude=lats,method='nearest')[month].plot(ax=ax3,x='latitude',y='depth',vmin=smin,vmax=smax,cmap=ocm.haline)
            ax3.set_title('SOSE')

            temp_plt = ax2.pcolormesh(x,y,rtemp_val[month],vmin=tmin,vmax=tmax,cmap=ocm.thermal)
            plt.colorbar(temp_plt,ax=ax2,label='temperature')
            ax2.set_title('ROMS')
            ax2.set_xlabel('latitude')
            ax2.set_ylabel('depth')

            salt_plt = ax4.pcolormesh(x,y,rsalt_val[month],vmin=smin,vmax=smax,cmap=ocm.haline)
            plt.colorbar(salt_plt,ax=ax4,label='salinity')
            ax4.set_title('ROMS')
            ax4.set_xlabel('latitude')
            ax4.set_ylabel('depth')

            if max_depth != None:
                for ax in [ax1,ax2,ax3,ax4]:
                    ax.set_ylim([-max_depth,0])

            plt.show()

        interact(plot,month=(0,11))
        
    elif mean == True:
        
        fig,axes = plt.subplots(2,2,figsize=(15,12))
        (ax1,ax2,ax3,ax4) = axes.flatten()

        sds.temperature.sel(longitude=slon,latitude=lats,method='nearest').mean('time').plot(ax=ax1,x='latitude',y='depth',vmin=tmin,vmax=tmax,cmap=ocm.thermal)
        ax1.set_title('SOSE')
        sds.salinity.sel(longitude=slon,latitude=lats,method='nearest').mean('time').plot(ax=ax3,x='latitude',y='depth',vmin=smin,vmax=smax,cmap=ocm.haline)
        ax3.set_title('SOSE')

        temp_plt = ax2.pcolormesh(x,y,np.mean(rtemp_val,axis=0),vmin=tmin,vmax=tmax,cmap=ocm.thermal)
        plt.colorbar(temp_plt,ax=ax2,label='temperature')
        ax2.set_title('ROMS')
        ax2.set_xlabel('latitude')
        ax2.set_ylabel('depth')

        salt_plt = ax4.pcolormesh(x,y,np.mean(rsalt_val,axis=0),vmin=smin,vmax=smax,cmap=ocm.haline)
        plt.colorbar(salt_plt,ax=ax4,label='salinity')
        ax4.set_title('ROMS')
        ax4.set_xlabel('latitude')
        ax4.set_ylabel('depth')
        
        if max_depth != None:
            for ax in [ax1,ax2,ax3,ax4]:
                ax.set_ylim([-max_depth,0])

        plt.show()


### define the plotting function
def compare_onshore(rds,sds,grid_coord,max_depth=None,tmin=-3,tmax=1,smin=33.8,smax=34.8):            

    print('define slice')
    [xi_min,eta_min,xi_max,eta_max,loc] = grid_coord

    

    rds_tr = rds.isel(eta_rho=slice(eta_min,eta_max+1),xi_rho=slice(xi_min,xi_max+1))

    x = rds_tr.lat_rho.values.squeeze()
    y = rds_tr.depth[0].values.squeeze()
    
    lons = rds_tr.lon_rho.to_masked_array().squeeze()
    lons[lons<0]+=360.0
    lats = rds_tr.lat_rho.to_masked_array().squeeze()
    
    lon_da = xr.DataArray(lons,dims='latitude')
    lat_da = xr.DataArray(lats,dims='latitude')
    
    print('calculate mean')
    temp_val = rds_tr.temp.mean('ocean_time').to_masked_array().squeeze()
    salt_val = rds_tr.salt.mean('ocean_time').to_masked_array().squeeze()
    
    print('plot')
    plt.close()
    
    fig,axes = plt.subplots(2,2,figsize=(15,10))
    (ax1,ax2,ax3,ax4) = axes.flatten()
    
    temp_plt = ax1.pcolormesh(x,y,temp_val,vmin=tmin,vmax=tmax,cmap=ocm.thermal)
    plt.colorbar(temp_plt,ax=ax1,label='Temperature in deg C')
    ax1.set_title('ROMS')
    ax1.set_xlabel('latitude')
    ax1.set_ylabel('depth')

    sds.temperature.sel(longitude=lon_da,latitude=lat_da,method='nearest').mean('time').plot(ax=ax2,vmin=tmin,vmax=tmax,cmap=ocm.thermal)
    ax2.set_title('SOSE')

    salt_plt = ax3.pcolormesh(x,y,salt_val,vmin=smin,vmax=smax,cmap =ocm.haline)
    plt.colorbar(salt_plt,ax=ax3,label='Salinity in PSU')
    ax3.set_title('ROMS')
    ax3.set_xlabel('latitude')
    ax3.set_ylabel('depth')

    sds.salinity.sel(longitude=lon_da,latitude=lat_da,method='nearest').mean('time').plot(ax=ax4,vmin=smin,vmax=smax,cmap=ocm.haline)
    ax4.set_title('SOSE')

    if max_depth != None:
        for ax in [ax1,ax2,ax3,ax4]:
            ax.set_ylim([-max_depth,0])
    
    plt.show()        

### define the plotting routine
def plot_cavity(rds,grid_coord,max_depth=None,tmin=-3,tmax=1,smin=33.8,smax=34.8,mean=True):            

    print('define slice')
    [xi_min,eta_min,xi_max,eta_max,loc] = grid_coord

    plt.close()

    rds_tr = rds.isel(eta_rho=slice(eta_min,eta_max+1),xi_rho=slice(xi_min,xi_max+1))
    
    x = rds_tr.lat_rho.values.squeeze()
    y = rds_tr.depth[0].values.squeeze()
    
    y[np.isnan(y)]=0.0
    
    if mean == True:
        
        fig,(ax1,ax2) = plt.subplots(nrows=2,figsize=(15,10))
        
        print('calculate mean')
        temp_val = rds_tr.temp.mean('ocean_time').to_masked_array().squeeze()
        salt_val = rds_tr.salt.mean('ocean_time').to_masked_array().squeeze()

        print('plot')
        temp_plt = ax1.pcolormesh(x,y,temp_val,vmin=tmin,vmax=tmax,cmap=ocm.thermal)
        plt.colorbar(temp_plt,ax=ax1,label='deg C')
        ax1.set_title('Temperature')
        ax1.set_ylabel('depth')

        salt_plt = ax2.pcolormesh(x,y,salt_val,vmin=smin,vmax=smax,cmap =ocm.haline)
        plt.colorbar(salt_plt,ax=ax2,label='PSU')
        ax2.set_title('Salinity')
        ax2.set_xlabel('latitude')
        ax2.set_ylabel('depth')

        if max_depth != None:
            for ax in [ax1,ax2]:
                ax.set_ylim([-max_depth,0])
        
        plt.show()
        
    elif mean == False:
        
        print('set up widget')
        def widget(month):
            
            fig,(ax1,ax2) = plt.subplots(nrows=2,figsize=(15,10))

            temp_val = rds_tr.temp.to_masked_array().squeeze()
            salt_val = rds_tr.salt.to_masked_array().squeeze()

            temp_plt = ax1.pcolormesh(x,y,temp_val[month],vmin=tmin,vmax=tmax,cmap=ocm.thermal)
            plt.colorbar(temp_plt,ax=ax1,label='deg C')
            ax1.set_title('Temperature')
            ax1.set_ylabel('depth')

            salt_plt = ax2.pcolormesh(x,y,salt_val[month],vmin=smin,vmax=smax,cmap =ocm.haline)
            plt.colorbar(salt_plt,ax=ax2,label='PSU')
            ax2.set_title('Salinity')
            ax2.set_xlabel('latitude')
            ax2.set_ylabel('depth')
            
            if max_depth != None:
                for ax in [ax1,ax2]:
                    ax.set_ylim([-max_depth,0])
            
            plt.show()
            
        interact(widget,month=(0,11))
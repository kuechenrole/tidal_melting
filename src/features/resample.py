from scipy import interpolate
import pyresample
import xarray as xr
import numpy as np

def resample(lon_s,lat_s,lon_t,lat_t,values):
    
    orig_def = pyresample.geometry.SwathDefinition(lons=lon_s,lats=lat_s)
    targ_def = pyresample.geometry.SwathDefinition(lons=lon_t,lats=lat_t)
    
    values_resampled = pyresample.kd_tree.resample_nearest(orig_def, values, \
        targ_def, radius_of_influence=500000, fill_value=None,nprocs=12)
    
    return values_resampled

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
    
    
    # Fill the mask of low resolution data with nearest neibghours and fill in known values on high res grid.
    if dim == 2:
        
        print('Fill in the mask of lr data')
        data = lr_da.values

        valid_mask = ~np.isnan(data)
        coords = np.array(np.nonzero(valid_mask)).T
        values = data[valid_mask]

        it = interpolate.NearestNDInterpolator(coords,values)

        filled = it(list(np.ndindex(data.shape))).reshape(data.shape)
        
        print('Resample to high resolution grid')
        hr_da[:,:] = resample(lr_grd['lon_'+gt].values,lr_grd['lat_'+gt].values,hr_grd['lon_'+gt].values,hr_grd['lat_'+gt].values,filled)
        
        # Fill with zeros where mask is present
        print('fill hr mask areas with fill value: ',fill_value)
        hr_da.values[hr_grd['mask_'+gt].values == 0] = fill_value
            
    if dim == 3:
        
        print('Fill in the mask of lr data, resample to high resolution grid and fill hr mask with fill value: ',fill_value)
        for k in np.arange(N):
            
            print('processing depth level: ',k)
            data = lr_da[k].values

            valid_mask = ~np.isnan(data)
            coords = np.array(np.nonzero(valid_mask)).T
            values = data[valid_mask]

            it = interpolate.NearestNDInterpolator(coords,values)

            filled = it(list(np.ndindex(data.shape))).reshape(data.shape)
    
            # Fill in known values on high res grid
            hr_da[k] = resample(lr_grd['lon_'+gt].values,lr_grd['lat_'+gt].values,hr_grd['lon_'+gt].values,hr_grd['lat_'+gt].values,filled)
            
            # Fill with zeros where mask is present
            print('fill hr mask areas with fill value: ',fill_value)
            hr_da[k].values[hr_grd['mask_'+gt].values == 0] = fill_value
            
    return hr_da

def low_to_high_frc(lr_da,lr_grd,hr_grd,gt,time_coord,fill_value=0.0):
    
    print('set up empty hr data array')
        


    N = time_coord.size
    dummy = np.tile(np.zeros(hr_grd['lon_'+gt].shape),(N,1,1))
    x = hr_grd['xi_'+gt]
    y = hr_grd['eta_'+gt]
    z = time_coord
    hr_da = xr.DataArray(dummy,coords=[z,y,x],dims=[time_coord.name,'eta_'+gt,'xi_'+gt])
    
    print('Fill in the mask of lr data, resample to high resolution grid and fill hr mask with fill value: ',fill_value)
    for k in np.arange(N):

        print('processing time: ',k)
        data = lr_da[k].values

        valid_mask = ~np.isnan(data)
        coords = np.array(np.nonzero(valid_mask)).T
        values = data[valid_mask]

        it = interpolate.NearestNDInterpolator(coords,values)

        filled = it(list(np.ndindex(data.shape))).reshape(data.shape)

        # Fill in known values on high res grid
        hr_da[k] = resample(lr_grd['lon_'+gt].values,lr_grd['lat_'+gt].values,hr_grd['lon_'+gt].values,hr_grd['lat_'+gt].values,filled)

        # Fill with zeros where mask is present

        hr_da[k].values[hr_grd['mask_'+gt].values == 0] = fill_value
            
    return hr_da
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title =([RunName,' surface heat/salt fluxes and wind stress']); % Make title from model name variable

disp('load air sea flux file (big file so load it once and subset later)')
load([RunName,'_air_sea_fluxes_daily.mat'])
disp('load uvstress')
load([RunName,'_ustress_daily.mat'])
load([RunName,'_vstress_daily.mat'])

%% 

disp(' ')
disp([' Creating the file : ',frcname])
disp(' ')


time_big = 365*(MaxYear-MinYear+1); %days
cycle = time_big;
Lpinfo=ncinfo(grdname,'lon_rho'); Lp = Lpinfo.Size(1);
Mpinfo=ncinfo(grdname,'lat_rho'); Mp = Mpinfo.Size(2);
L=Lp-1;
M=Mp-1;
%% 

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('clobber'));
cmode = bitor(cmode,netcdf.getConstant('Classic_Model'));

id = netcdf.create(frcname, cmode);
%
%  Create dimensions
%
    disp('Create dimension')

xi_u_dim = netcdf.defDim(id, 'xi_u', L);
eta_u_dim = netcdf.defDim(id, 'eta_u', Mp);
xi_v_dim = netcdf.defDim(id, 'xi_v', Lp);
eta_v_dim = netcdf.defDim(id, 'eta_v', M);
xi_rho_dim = netcdf.defDim(id, 'xi_rho', Lp);
eta_rho_dim = netcdf.defDim(id, 'eta_rho', Mp);
xi_psi_dim = netcdf.defDim(id, 'xi_psi', L);
eta_psi_dim = netcdf.defDim(id, 'eta_psi', M);
sms_time_dim = netcdf.defDim(id, 'sms_time', time_big);
shf_time_dim = netcdf.defDim(id, 'shf_time', time_big);
swf_time_dim = netcdf.defDim(id, 'swf_time', time_big);
%
%  Create variables and attributes
%
    disp('Create variable attribute')
sms_time_id = netcdf.defVar(id, 'sms_time', 'double', sms_time_dim);
netcdf.putAtt(id, sms_time_id, 'long_name', 'surface momentum stress time');
netcdf.putAtt(id, sms_time_id, 'units', 'day');
netcdf.putAtt(id, sms_time_id, 'cycle_length', cycle);

shf_time_id = netcdf.defVar(id, 'shf_time', 'double', shf_time_dim);
netcdf.putAtt(id, shf_time_id, 'long_name', 'surface heat flux time');
netcdf.putAtt(id, shf_time_id, 'units', 'day');
netcdf.putAtt(id, shf_time_id, 'cycle_length', cycle);

swf_time_id = netcdf.defVar(id, 'swf_time', 'double', swf_time_dim);
netcdf.putAtt(id, swf_time_id, 'long_name', 'surface freshwater flux time');
netcdf.putAtt(id, swf_time_id, 'units', 'day');
netcdf.putAtt(id, swf_time_id, 'cycle_length', cycle);

sustr_id = netcdf.defVar(id, 'sustr', 'double', [xi_u_dim, eta_u_dim sms_time_dim]);
netcdf.putAtt(id, sustr_id, 'long_name', 'surface u-momentum stress');
netcdf.putAtt(id, sustr_id, 'units', 'Newton meter-2');
netcdf.putAtt(id, sustr_id, 'missing_value', -1e34);
netcdf.putAtt(id, sustr_id, '_FillValue', -1e34);

svstr_id = netcdf.defVar(id, 'svstr', 'double', [xi_v_dim, eta_v_dim sms_time_dim]);
netcdf.putAtt(id, svstr_id, 'long_name', 'surface v-momentum stress');
netcdf.putAtt(id, svstr_id, 'units', 'Newton meter-2');
netcdf.putAtt(id, svstr_id, 'missing_value', -1e34);
netcdf.putAtt(id, svstr_id, '_FillValue', -1e34);

shflux_id = netcdf.defVar(id, 'shflux', 'double', [xi_rho_dim, eta_rho_dim shf_time_dim]);
netcdf.putAtt(id, shflux_id, 'long_name', 'surface net heat flux');
netcdf.putAtt(id, shflux_id, 'units', 'Watts meter-2');
netcdf.putAtt(id, shflux_id, 'missing_value', -1e34);
netcdf.putAtt(id, shflux_id, '_FillValue', -1e34);

swflux_id = netcdf.defVar(id, 'swflux', 'double', [xi_rho_dim, eta_rho_dim swf_time_dim]);
netcdf.putAtt(id, swflux_id, 'long_name', 'surface freshwater flux (E-P)');
netcdf.putAtt(id, swflux_id, 'units', 'centimetre day-1');
netcdf.putAtt(id, swflux_id, 'missing_value', -1e34);
netcdf.putAtt(id, swflux_id, '_FillValue', -1e34);
netcdf.putAtt(id, swflux_id, 'positive', 'net evaporation');
netcdf.putAtt(id, swflux_id, 'negative', 'net precipitation');

%
% Create global attributes
%
    disp('Create global attribute')
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', title);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'clim_file', frcname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'ROMS forcing file');
netcdf.endDef(id);
%% 


    disp('Write time variable')

swft = [1:cycle]; %data at each day of cycle
shft = [1:cycle];

netcdf.putVar(id, sms_time_id, swft);
netcdf.putVar(id, shf_time_id, shft);
netcdf.putVar(id, swf_time_id, swft);

    disp('write salt and heat flux variables')
shf = shfluxGrid; clear shfluxGrid
swf = ssfluxGrid; clear ssfluxGrid
su  = u_stress; clear u_stress
sv  = v_stress; clear v_stress

    disp('inpainting nans in wind fields...')
%addpath('/ds/projects/iomp/matlab_scripts')
for j = 1:time_big;
  shf(j,:,:)= inpaint_nans(squeeze( shf(j,:,:)),1);
  swf(j,:,:)= inpaint_nans(squeeze( swf(j,:,:)),1);
  su(j,:,:) = inpaint_nans(squeeze(su(j,:,:)),1);
  sv(j,:,:) = inpaint_nans(squeeze(sv(j,:,:)),1);
if ~rem(j,time_big/20)
disp([num2str(j/time_big*100),'%'])
end
end
%rmpath('/ds/projects/iomp/matlab_scripts')


    disp('correct for large summer heat flux')

% Correct Heat flux for large summer values (as per NCEP2 heat flux data)
shf(shf > 0) = shf(shf > 0)*0.5;
%swf(swf <=0) = -1e-6;

refSalt = 34.4; % Reference salinity


netcdf.putVar(id, shflux_id, permute(shf,[3 2 1]));% (W/m^2)
netcdf.putVar(id, swflux_id, permute(-swf./refSalt*100,[3 2 1]));% kg/m^2->cm day^-1
%netcdf.putVar(id, ssflux_id, permute(-swf./refSalt/(60*60*24),[3 2 1])); %kg/m^2 -> m/s
    disp('Write wind variables')

netcdf.putVar(id, sustr_id, permute((su(:,:,1:end-1)+su(:,:,2:end))./2,[3 2 1]));
netcdf.putVar(id, svstr_id, permute((sv(:,1:end-1,:)+sv(:,2:end,:))./2,[3 2 1]));

netcdf.close(id);

   disp('DONE')
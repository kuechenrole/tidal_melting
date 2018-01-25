%% ERA_interim_misom_grid_stress_annual.m
%
% Load WIND speed/dir from ERA-interim, and convert to stress at the surface. Specifically, take ERA-interim daily data, and output as daily stress.
% E Cougnon - Wrote original script for loading ERA data (July 2014)
% D Gwyther - adapted for 1-daily forcing, some corrections and alterations to fit in my BC creator framework (dec 2014)
% D Gwyther - Update to remove old dependencies (apr 2016)
%%

%addpath('/ds/projects/iomp/matlab_scripts')
% Load model grid (grdname from make_sbc.m):
lat_rho=ncread(grdname,'lat_rho')';
lon_rho=ncread(grdname,'lon_rho')';
mask_rho=ncread(grdname,'mask_rho')';
angle=ncread(grdname,'angle')';

NumYears = MaxYear-MinYear+1;
SamRate = 2;% sample rate of data (default for ERA_interim is 1/day)

Imin_ERAi = Imin_wind; %lonmin
Imax_ERAi = Imax_wind;
Jmin_ERAi = Jmin_wind; %latmin
Jmax_ERAi = Jmax_wind;

u10Path = [external_dir,'/era_interim/ERA_Interim_1992_2011.2daily.uwinds.nc'] %'/home/ubuntu/katabatic/obs/ERA_Interim/ERA_Interim_1992_2011.2daily.uwinds.nc'
v10Path = [external_dir,'/era_interim/ERA_Interim_1992_2011.2daily.vwinds.nc'] %'/home/ubuntu/katabatic/obs/ERA_Interim/ERA_Interim_1992_2011.2daily.vwinds.nc'
disp('loading ERA-interim data')
u10=ncread(u10Path,'u10',[Imin_ERAi Jmin_ERAi 1],[Imax_ERAi-Imin_ERAi+1 Jmax_ERAi-Jmin_ERAi+1 Inf]);
u10 = permute(u10,[3 2 1]);
v10=ncread(v10Path,'v10',[Imin_ERAi Jmin_ERAi 1],[Imax_ERAi-Imin_ERAi+1 Jmax_ERAi-Jmin_ERAi+1 Inf]);
v10 = permute(v10,[3 2 1]);
longitude=double(ncread(u10Path,'longitude',[Imin_ERAi],[Imax_ERAi-Imin_ERAi+1]));
for lonInd=1:size(longitude,1);
    if longitude(lonInd) > 180.0
     longitude(lonInd) = longitude(lonInd)-360.0;
    end
end
latitude=double(ncread(v10Path,'latitude', [Jmin_ERAi],[Jmax_ERAi-Jmin_ERAi+1]));
uwndall=[];
vwndall=[];
u_stress_All = [];
v_stress_All = [];
u_index=1;
v_index=1;
%% 

LeapYears = [1992:4:2040]; %leap years til 2040
for YearInd = MinYear:MaxYear; 
clear uwnd vwnd uw_stress vw_stress signu signv
 if any(YearInd == LeapYears)
          % allocate space for matrix for daily data for one year
  uwnd = nan(366*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  vwnd = nan(366*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  uw_stress = nan(366,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  vw_stress = nan(366,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);          
 else
          % allocate space for matrix for daily data for one year
  uwnd = nan(365*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  vwnd = nan(365*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  uw_stress = nan(365,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  vw_stress = nan(365,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
 end

 uwnd(:,:,:) = squeeze(u10(u_index:u_index+size(uwnd,1)-1,:,:));
 vwnd(:,:,:) = squeeze(v10(v_index:v_index+size(vwnd,1)-1,:,:));

   signu = sign(uwnd);
   signv = sign(vwnd);

   rhoAir = 1.3;
   Cd = 1.4e-3;

   taux = rhoAir*Cd.*uwnd.^2.*signu;
   tauy = rhoAir*Cd.*vwnd.^2.*signv;

   % Want daily data, so far data are 2 daily -- sample one every 2 data to get the data
   % taken from 12.00 or average 00:00 and 12:00
if SamRate == 1 %daily data at 12:00 only
       uw_stress=taux;
       vw_stress=tauy;
elseif SamRate == 2 % daily data at 00:00 and 12:00
   k=1;
   for i=1:2:length(uwnd)-1;
       uw_stress(k,:,:) = 0.5*(taux(i,:,:)+taux(i+1,:,:));
       vw_stress(k,:,:) = 0.5*(tauy(i,:,:)+tauy(i+1,:,:));	
       k=k+1;
   end
end
%% 

   AISuw_stress=[];
   AISvw_stress=[];
   

deg2rad = pi/180.0
roms_x = -(lat_rho+90).*cos(lon_rho*deg2rad+pi/2);
roms_y = (lat_rho+90).*sin(lon_rho*deg2rad+pi/2);
[lon_era,lat_era] = meshgrid(longitude,latitude);
era_x = -(lat_era+90).*cos(lon_era*deg2rad+pi/2);
era_y = (lat_era+90).*sin(lon_era*deg2rad+pi/2);


   % Interpolate each daily data to ROMS grid
disp(['Interpolating and rotating daily data to ROMS grid for ',num2str(YearInd)])
   for ii = 1:size(uw_stress,1);
       disp(['processing day: ',num2str(ii)]) 
       ustress_nonrot = griddata(era_x,era_y,squeeze(uw_stress(ii,:,:)),roms_x,roms_y,'linear');
       vstress_nonrot = griddata(era_x,era_y,squeeze(vw_stress(ii,:,:)),roms_x,roms_y,'linear');
       [urot,vrot]=rotate_vec(ustress_nonrot,vstress_nonrot,angle,0);
       AISuw_stress(ii,:,:) = urot;
       AISvw_stress(ii,:,:) = vrot;
   end    
    
   disp('resampling done')
   u_stress=AISuw_stress;
   v_stress=AISvw_stress;
   %u_stress_All(u_index:u_index+size(uw_stress)-1,:,:) = AISuw_stress;
   %v_stress_All(v_index:v_index+size(uw_stress)-1,:,:) = AISvw_stress;


   disp(['Saving u_stress and v_stress for ' num2str(YearInd) ''])
   % save u and v component with the model grid, and a vector rotation
   nameval=[RunName,'_ustress_daily.mat'];
   save(nameval,'u_stress','-v7.3'); 
   nameval=[RunName,'_vstress_daily.mat'];
   save(nameval,'v_stress','-v7.3');
   disp(['' num2str(YearInd) ' Saved.']),

   %u_index = u_index+size(uw_stress,1);
   %v_index = v_index+size(uw_stress,1);
end

% disp(['Saving u_stress and v_stress for ' num2str(MinYear) ' to ' num2str(MaxYear) ' with backup leap year data'])
%   % save u and v component with the model grid, and a vector rotation
%   nameval=['ustress_grid_model_',num2str(MinYear),'_',num2str(MaxYear),'withleapyear.mat'];
%   save(nameval,'u_stress_All','-v7.3');
%   nameval=['vstress_grid_model_' num2str(MinYear),'_',num2str(MaxYear),'withleapyear.mat'];
%   save(nameval,'v_stress_All','-v7.3');


%LeapYear = [1992:4:2040]; %leap years til 2040
%disp('Removing leap year data: cut and remove data on feb-29 leap years')
% METHOD 1
%ly=zeros([1 366]); ly(60)=1;
%nly=zeros([1 365]);
%ly_index=ismember([MinYear:MaxYear],LeapYear); 
%feb29_index=[];
%for ii=1:length(ly_index)
%if ly_index(ii)==1
%feb29_index=[feb29_index,ly];
%elseif ly_index(ii)==0
%feb29_index=[feb29_index,nly];
%end
%end
% METHOD 2
%shfluxClima_tmp = nans(366,MaxYear-MinYear+1,x,y);
%DayLoop=1;
%for yy = 1:MaxYear-MinYear+1
%for dd = 1:365
%if any(yy+MinYear-1==LeapYear) & dd==60;
%shfluxClima_tmp(dd,yy,x,y) = shfluxGrid(DayLoop+1,:,:);
%shfluxClima_tmp(366,yy,x,y) = shfluxGrid(DayLoop,:,:);
%DayLoop=DayLoop+1;
%else 
%shfluxClima_tmp(dd,yy,x,y) = shfluxGrid(DayLoop,:,:);
%end
%DayLoop=DayLoop+1;
%end
%end
%shfluxClima=nanmean(shfluxClima_tmp,1);
% METHOD 3
%feb29_index=zeros([1 size(u_stress_All,1)]);
%FLYpos = find(ismember([MinYear:MaxYear],LeapYear),1);
%for ii=0:length(find(ismember([MinYear:MaxYear],LeapYear)))-1
%feb29_index( ((FLYpos-1)*365+60) + ii*(306+(365*3)+60))=1;
%end
%uwndClima_tmp=u_stress_All; clear u_stress_All
%vwndClima_tmp=v_stress_All; clear v_stress_All
%for ii=find(feb29_index)'%loop through feb-29 indices
%uwndClima_tmp(ii,:,:)=[]; %remove feb29 values
%vwndClima_tmp(ii,:,:)=[]; %remove feb29 values
%end

%u_stress_All = uwndClima_tmp;
%v_stress_All = vwndClima_tmp;

% disp(['Saving u_stress and v_stress for ' num2str(MinYear) ' to ' num2str(MaxYear)])
%   % save u and v component with the model grid, and a vector rotation
%   nameval=['ustress_grid_model.mat'];
%   save(nameval,'u_stress_All','-v7.3');
%   nameval=['vstress_grid_model.mat'];
%   save(nameval,'v_stress_All','-v7.3');


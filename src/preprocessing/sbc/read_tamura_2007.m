%% READ_TAMURA_DAILY
% Read and save daily ERA-interim data from T. Tamura and save to mat file for ROMS surface forcing file creation
%
% Created by B K. Galton-Fenzi 
% Adapted to Daily forcing by E. Cougnon at /ds/projects/iomp/obs/Tamura_air_sea_fluxes/daily/ (May 2013)
% Poked around a bit by D Gwyther at /ds/projects/iomp/totten/ana/dgwyther/netcdf/forcing_creator/sbc/read_tamura_daily.m (March and Dec 2014)
% major changes by D Gwyther, Jan 2016
% merged to new file for just 1995 (defined as a normal clima year) by D Gwyther, Jul 2016
%%%%%%

lon_rho=ncread(grdname,'lon_rho')';
lat_rho=ncread(grdname,'lat_rho')';

shfluxtmp = nan(xmax-xmin+1,ymax-ymin+1,365); 
ssfluxtmp = nan(xmax-xmin+1,ymax-ymin+1,365); 
shflux = nan(365,xmax-xmin+1,ymax-ymin+1);
ssflux = nan(365,xmax-xmin+1,ymax-ymin+1);

%% Read in land mask
fid3=fopen([external_dir,'/tamura/EASE_landmask_H.data'],'r');
landmaskNaN = reshape(fread(fid3,721*721*1,'float32=>double'),721,721);
landmaskNaN(landmaskNaN==0)=NaN;

ii = 1;
Month = ['jan';'feb';'mar';'apr';'may';'jun';'jul';'aug';'sep';'oct';'nov';'dec'];
DaysPerMonth=[31,28,31,30,31,30,31,31,30,31,30,31];
%% Read in data:
for YearInd = 2007
 for MonthInd = 1:12; %Keep track of current month
    dms = DaysPerMonth(MonthInd);

    display(['Reading in GrADS data ' num2str(YearInd),' ',Month(MonthInd,:),' ...']);
    fid = fopen([external_dir,'/tamura/TSDM2hb_' num2str(YearInd),'_',Month(MonthInd,:),'.data']);
    
    tmp = reshape(fread(fid,721*721*dms*6,'float32=>double'),721,721,6,dms); 
    tmp = bsxfun(@times,tmp,landmaskNaN); % land mask
    
    
    shfluxtmp(:,:,ii:ii+dms-1) = squeeze(tmp(xmin:xmax,ymin:ymax,1,:)); %[W/m^2] 4:NCEP 1:ERA-interim
    ssfluxtmp(:,:,ii:ii+dms-1) = squeeze(tmp(xmin:xmax,ymin:ymax,3,:)); %[Kg/m^2] 6:NCEP 3:ERA-interim
    
  ii = ii+dms;
 end
save YearInd.txt YearInd -ascii
end

DayNumber = ii-1;% save the day index to a more descriptive Varname (don't forget to remove 1 to account forfirst day being indexed)

if 1 % do inpainting
 for iii = 1:size(shfluxtmp,3);
    disp(['interpolating NaNs in fluxes...',num2str(iii/size(shfluxtmp,3)*100),'% done'])
    shflux(iii,:,:) = inpaint_nans(squeeze(shfluxtmp(:,:,iii)),2); % Heat: Watts/m^2
    ssflux(iii,:,:) = inpaint_nans(squeeze(ssfluxtmp(:,:,iii)),2); % (m/month)/(s/month)*(kg/m^2)  cm/day
 end
elseif 0
shflux = permute(shfluxtmp,[3 1 2]);
ssflux = permute(ssfluxtmp,[3 1 2]);
end

%%

disp('Saving Heat/Salt fluxes from Takeshi''s data')
save([RunName,'_Takeshi_subset_daily.mat'],'shfluxtmp','ssfluxtmp')
disp('Saved, now gridding to model grid')

fid2 = fopen([external_dir,'/tamura/latlon.data']);
B = fread(fid2,721*721*2,'float32=>double');
BB = reshape(B,721,721,2);

display('Making tidy for ROMS grid...');
AllLon = squeeze(BB(xmin:xmax,ymin:ymax,2));
AllLat = squeeze(BB(xmin:xmax,ymin:ymax,1));

ssfluxGrid = nan(DayNumber,size(lat_rho,1),size(lat_rho,2));
shfluxGrid = nan(DayNumber,size(lat_rho,1),size(lat_rho,2));

deg2rad = pi/180.0;
roms_x = -(lat_rho+90).*cos(lon_rho*deg2rad+pi/2);
roms_y = (lat_rho+90).*sin(lon_rho*deg2rad+pi/2);
tamura_x = -(AllLat+90).*cos(AllLon*deg2rad+pi/2);
tamura_y = (AllLat+90).*sin(AllLon*deg2rad+pi/2);

for j = 1:DayNumber;
    ssfluxGrid(j,:,:) = griddata(tamura_x,tamura_y,squeeze(ssflux(j,:,:)),roms_x,roms_y,'nearest');
    shfluxGrid(j,:,:) = griddata(tamura_x,tamura_y,squeeze(shflux(j,:,:)),roms_x,roms_y,'nearest');
    if ~rem(j,round(DayNumber/10)), disp(['gridding ',num2str(j/DayNumber*100),' done.']), end
end


disp('Saving gridded Heat/Salt fluxes for 2007')
save([RunName,'_air_sea_fluxes_daily.mat'],'shfluxGrid','ssfluxGrid','-v7.3')
disp('Saved.')



smooth = 0;
PlotFigs = 0;
%addpath(genpath('/home/ubuntu/iceOceanVolume/matlab_tools'))
%addpath(genpath('/home/ubuntu/iceOceanVolume/WholeAntarcticModel/grid_generator/'))
%addpath(genpath('/home/ubuntu/iceOceanVolume/WholeAntarcticModel/Preprocessing_tools/'))
%addpath(genpath('/home/ubuntu/iceOceanVolume/ant_mf/bedmap2'))
%addpath(genpath('/home/ubuntu/iceOceanVolume/Code/roms_matlab'))
proj_dir = getenv('proj_dir');
extern_dir = [proj_dir,'/data/preprocessing/external'];
addpath(genpath('../matlab_tools/'))
addpath(genpath('../../../data/preprocessing/external/bedmap2'))

load hot_cold_white
%% Model domain at mesh resolution (mr) in km:
mr =  5

%establish domain size of roms mesh (South pole is at km 0,0):
[Cx Cy] = meshgrid([-4500:mr:4500],[-4500:mr:4500]);
% set up arrays for grid ice mask (Cim), bed height (Cb), ice thickness (Ci),
% water/land mask (Cwm), latitude (Cla), longitude (Clo) 
Cim = zeros(size(Cx));
Cb  = nan(size(Cx));
Ci  = nan(size(Cx));
Cwm = ones(size(Cx));

%establish longitude and latitude from polar stereo projection:
Cla = nan(size(Cx));
Clo = nan(size(Cx));
for i = 1:size(Cx,1);
    [t1, t2] = inverse_polar_stereo(Cx(i,:),Cy(i,:),0,0,1,-71);
    Cla(i,:) = t1';
    Clo(i,:) = t2';
end

%establish grid spacing variables: 
dx = mr.*1000.*ones(size(Cx));
dy = dx;
%%
load bedmap2

%rename bedmap x,y and bed elevation to bm_* for consistency
bm_x=x;
bm_y=y;
bm_b=-elev_bed;

%calculate bedmap latitude and longitude arrays:
%bm_lo = nan(6667,6667); bm_la = bm_lo;
%for i = 1:6667;
%[t1 t2] = inverse_polar_stereo(bm_x(i,:),bm_y(i,:),0,0,1,-71);
%bm_la(i,:) = t1';
%bm_lo(i,:) = t2';
%end

%modify Bedmap (B) masks to roms format:
%res = 10;
% ice draft roms format is negative elevation and zeros instead of Nans,
% where no ice
bm_draft = -(thick-elev_surf);
bm_draft(isnan(bm_draft) == 1) = 0;
% ocean mask gets 1 where open ocean (ice shelf excluded) and zero where land/ice shelf for now 
bm_wm = thick;
bm_wm(isnan(thick) == 1) = 1;
bm_wm(isnan(thick) == 0) = 0;
bm_wm(find(double(mask_rock) == 1)) = 0;
% ice mask (already just ice shelf == True) and mask of lake vostock
bm_im = double(mask_ice);
mask_vostok = double(mask_vostok);


%% Add in mask_vostok

% find indices (ii,jj), where left bottom corner of lake vostock mask fits into bigger ice mask 
% in x-direction
ii = find(bm_x(1,:) - x_mask_vostok(1,1) == 0)
%in y-direction
jj = find(bm_y(:,1) - y_mask_vostok(1,end) == 0)
% update values of ice mask from slected start indices to end of vostock
% mask
bm_im(jj:jj+size(mask_vostok,1)-1, ii:ii+size(mask_vostok,2)-1) = mask_vostok;
% set ice mask to zero, where bedmap has rock outcrops (rock exposure to
% atmosphere?)
bm_im(find(double(mask_rock) == 1)) = 0;
% add ice shelf as ones to ocean mask (roms rho-points)
bm_wm(find(double(bm_im) == 1)) = 1;
%% merge bedmap variables roms grid (two options):

% first find indices (ii,jj), where left bottom corner of grid domain fits into bigger bedmap domain 
% this is done by comparing the distance-from-pole arrays of bedmap and grid
%ii = find(x(1,:) - Cx(1,1) == 0)
%jj = find(y(:,1) - Cy(1,end) == 0)
% subset ice mask from selected start indices to end of grid size in steps of mr resolution
% CAUTION: this approach assumes that grid domain extend is in all
% directions smaller than bedmap domain extend.
%Cim = Cim(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
%Cb = -elev_bed(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
%Ci = draft(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
%Cwm =  Cwm(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
%Cla =  bm_la(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
%Clo =  bm_lo(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);

% first, find indices (ii,jj), where left bottom corner of grid domain fits into bigger bedmap domain 
% this is done by comparing the distance-from-pole arrays of bedmap and grid
ii = find(bm_x(1,4) - Cx(1,:) == 0)
jj = find(bm_y(4,1) - Cy(:,1) == 0);
% merge bedmap variables (small grid) in predefined roms grid (large grid):
Cim(jj:jj+size(bm_x(4:mr:end,4:mr:end),1)-1,ii:ii+size(bm_x(4:mr:end,4:mr:end),2)-1) = bm_im(4:mr:end,4:mr:end);
%Cb = -elev_bed(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
Cb(jj:jj+size(bm_x(4:mr:end,4:mr:end),1)-1,ii:ii+size(bm_x(4:mr:end,4:mr:end),2)-1) = bm_b(4:mr:end,4:mr:end);
%Cb(isnan(Cb) == 1) = 5000;
%Ci = draft(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
Ci(jj:jj+size(bm_x(4:mr:end,4:mr:end),1)-1,ii:ii+size(bm_x(4:mr:end,4:mr:end),2)-1) = bm_draft(4:mr:end,4:mr:end);
%Cwm =  Cwm(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
Cwm(jj:jj+size(bm_x(4:mr:end,4:mr:end),1)-1,ii:ii+size(bm_x(4:mr:end,4:mr:end),2)-1) = bm_wm(4:mr:end,4:mr:end);
%Cla =  bm_la(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
%Cla(jj:jj+size(bm_x(4:mr:end,4:mr:end),1)-1,ii:ii+size(bm_x(4:mr:end,4:mr:end),2)-1) = bm_la(4:mr:end,4:mr:end);
%Clo =  bm_lo(jj:mr:jj+(size(Cx,1)*mr)-1,ii:mr:ii+(size(Cy,2)*mr)-1);
%Clo(jj:jj+size(bm_x(4:mr:end,4:mr:end),1)-1,ii:ii+size(bm_x(4:mr:end,4:mr:end),2)-1) = bm_lo(4:mr:end,4:mr:end);
%define water column thickness as bed elevation + (negative ice draft)
Cwct = Cb+Ci;
% put NaNs back into ice and water mask
CimNaN = Cim;
CimNaN(Cim == 0) = NaN;
CwmNaN = Cwm;
CwmNaN(Cwm == 0) = NaN;

DTOR=pi/180.0;
angle=Clo.*DTOR;

% combine variables in matlab struct s 
s.x = Cx;
s.y = Cy;
s.dx = dx;
s.dy = dy;
s.lon = Clo;
s.lat = Cla;
s.mw = Cwm;
s.mi = Cim;
s.zb = Cb;
s.zd = Ci;
s.angle = angle;%zeros(size(Cy)); % angle should just be just angle = - longitude*pi/180
% defining minimum and maximum depth of roms grid
s.clipping_depths(1) = 20;
s.clipping_depths(2) = 6000;

%% Check point plots
if PlotFigs == 1;
h1 = figure('visible','off'),
figure,
flat_pcolor(s.zd),
%flat_pcolor(x(1:res:end,1:res:end),y(1:res:end,1:res:end),elev_bed(1:res:end,1:res:end));shading flat,
colormap(hot_cold_white),
%hold on,
%resg2 = 10;
%plot(Cx(:,1:resg2:end),Cy(:,1:resg2:end),'k');
%hold on,
%plot(Cx(1:resg2:end,:)',Cy(1:resg2:end,:)','k');
colorbar();
axis square
%caxis([-1500 1500])
%saveas(h1,'grd_20170809.png','png');
%!cp -p grd_20150513.png ~/public_html/WholeAntarcticModel/.
end


%% Fill up far north regions, with RTOPO (Smith and Sandwell bathymetry)
% Use rtopo bathymetry south of 30S (subset of global 30s bedrock variable)
RTOPO = '../../../data/preprocessing/external/rtopo/RTopo-2.0.1_30sec_bedrock_topography_S30.nc'
%RTOPO = [extern_dir,'/rtopo/RTopo-2.0.1_30sec_bedrock_topography_S30.nc'];
% load lon lat and bathy from RTopo and generate lat-lon mesh for interpolation
lon_rtopo = ncread(RTOPO,'lon');
lat_rtopo = ncread(RTOPO,'lat');
bathy_rtopo = ncread(RTOPO,'bedrock_topography')';
[lon_rtopo lat_rtopo] = meshgrid(lon_rtopo,lat_rtopo);
%ii = find(s.zb(:) == s.zb(1,1));
% find indices of grid bathy where there is no topo yet
ii = isnan(s.zb);
% interpolate RTopo to grid bathy to these points (interpolation based on lat lon coordinates)
out = interp2(lon_rtopo,lat_rtopo,bathy_rtopo,s.lon(ii),s.lat(ii));
%update grid bathy at the interpolayted locations
s.zb(ii) = -out;
% Generate a smooth transition between Rtopo and Bedmap
% Choose suture zone thickness = 3*mr;
sr = round((mr/110*3)*10)/10;
% select points around the edge of Bedmap (thickness depends on suture zone thickness) 
ii = find(s.lat(:) > -60.05 & s.lat(:) < (-60.05+sr)); % Crop suture zone and smoothly fill
% set these points to Nan
s.zb(ii) = NaN;
% interpolate these points from known neighbours
out = inpaint_nans(s.zb);
% update grid bathy at suture zone
s.zb = out;
% update water mask to include far north land
ii = find(s.zb<0.0 & s.lat > -60.0);
s.mw(ii) = 0;

% update ice draft to 0.0 at far north regions
ii = find(s.lat > -60.0);
s.zd(ii) = 0;
%% Check point plots
if PlotFigs == 1
h2 = figure(),
%figure('units','normalized','outerposition',[0 0 1 1])
flat_pcolor(s.zb)
colormap(hot_cold_white)
colorbar()
hold on
contour(s.lat,[-60.05 -60.05],'LineWidth',0.1,'LineColor','k')
x_vals=[-4300:500:-3300,-3000,3300:500:4300];
y_vals = [-3700:500:-2700,2600:500:3600];
%vals = [-3600:mr*10:3600];
contour(s.x,x_vals,'ShowText','on','LineColor','k')
contour(s.y,y_vals,'ShowText','on','LineColor','k')
%plot(Cx(:,:),Cy(:,:),'k');
%hold on, plot(Cx(:,:)',Cy(:,:)','k'););
%caxis([-1500,1500]);20170809
%saveas(h2,'waom10_grd_domains.png','png');
axis square,
hold off;
end

%% Define different domain sizes
x_orig = [-3000,3300];
y_orig = [-2700,2600];

x_s = [-3300,3300];
y_s = [-2700,2600];

x_m = [-3800,3800];
y_m = [-3200,3100];

x_l = [-4300,4300];
y_l = [-3700,3600];

%subset the same way Ben did befor with find lower corner and subset up to
%upper corner!

i = find(s.x(1,:) == x_l(1))
ii = find(s.x(1,:) == x_l(2))

j = find(s.y(:,1) == y_l(1))
jj = find(s.y(:,1) == y_l(2))

s.x = s.x(j:jj,i:ii);
s.y = s.y(j:jj,i:ii);
s.dx = s.dx(j:jj,i:ii);
s.dy = s.dy(j:jj,i:ii);
s.lon = s.lon(j:jj,i:ii);
s.lat = s.lat(j:jj,i:ii);
s.mw = s.mw(j:jj,i:ii);
s.mi = s.mi(j:jj,i:ii);
s.zb = s.zb(j:jj,i:ii);
s.zd = s.zd(j:jj,i:ii);
s.angle = s.angle(j:jj,i:ii);

s.mwNaN = CwmNaN(j:jj,i:ii);
s.miNaN = CimNaN(j:jj,i:ii);

% check corners again + how many cells? 528 or 529?
%% CREATE ROMS GRID
projection = 'Polar stereographic, lat0 = -71 S';
% ROMS needs radians.
grid_x = s.x;
grid_y = s.y;
[m, n] = size(grid_x);
geogrid_lon = s.lon;
geogrid_lat = s.lat;
geometry{1} = s.dx;
geometry{2} = s.dy;
mask = s.mw(2:end,2:end);   % land = 1; water = 0.ii = find(x_s(1) <= Cx <= x_s(2));
imask = s.mi(2:end,2:end);
if ~isequal(size(mask), size(grid_x)-1)
if ~isempty(mask)
disp(' ## Wrong size mask.')
end
mask = zeros(m-1, n-1);
end
if ~isequal(size(imask), size(grid_x)-1)
if ~isempty(imask)
disp(' ## Wrong size mask.')
end
imask = zeros(m-1, n-1);
end
%mask = ~mask;
%land = mask;
%water = ~land;
%imask = imask;
projection = 'Stereographic';
ang = s.angle;   % ROMS needs radians.


min_depth = s.clipping_depths(1);
min_depth0 = 0.1; % 10 cm
max_depth = s.clipping_depths(2);


% Smoothing of bathy to reduce PGF errors:
CwmFALSE = ones(size(s.mw));
bt = s.zb;
dt = s.zd;
wct = bt+dt;


% Crop mountains
dt(bt < min_depth0) = 0;
bt(bt < min_depth0) = min_depth0;
wct = bt+dt;
dt(wct < 0) = -bt(wct < 0) + min_depth0;

% Preserve subglacial regions
dt(wct == 0) = -bt(wct == 0) + min_depth0;

wct = bt + dt;

rx1in = 0.3

if smooth == 1
    disp('smoothing activated');
    [wctout]=smooth_bath(wct,CwmFALSE,4,rx1in,150);
    disp("smooth wct ok");
    [bathymetry]=smooth_bath(bt,CwmFALSE,4,rx1in,150);
    disp("smooth bathy ok");
else
    disp('smoothing deactivated');
    wctout = wct;
    bathymetry = bt;
end

ice_draft = wctout - bathymetry;
wct = wctout;

s.mis = s.mw;
s.mis ~= s.mw;
s.mis(s.mi == 1) = 0;
s.mis = double(~s.mis);


ice_draft(s.mw == 0) = -bathymetry(s.mw == 0) + min_depth0;
wct = bathymetry + ice_draft;
ii = find(wct < min_depth0 & s.mw == 1);
bathymetry(ii) = -ice_draft(ii) + min_depth;

% Make Water column thickness satisfys min_depth
s.mw00 = s.mw;
s.mw00(s.mw == 0) = NaN;
% Check no positive ice draft values
ice_draft(ice_draft > 0) = 0;
ice_draft(s.mw == 0) = -bathymetry(s.mw == 0) + min_depth0;
wct = bathymetry + ice_draft;
disp("wct min thickness ok");


mindxy = mr*1000; % 1kmmin(dx,dy);
CFL = mindxy./abs(sqrt(9.81.*wct)).*s.mwNaN;
%
[out IndminCFL] = nanmin(CFL(:)),
disp("nan_min cfl ok");


theInterpFcn = 'interp2';
theInterpMethod = 'nearest';
grid_x = feval(theInterpFcn, grid_x, 1, theInterpMethod);
grid_y = feval(theInterpFcn, grid_y, 1, theInterpMethod);
geogrid_lon = feval(theInterpFcn, geogrid_lon, 1, theInterpMethod);
geogrid_lat = feval(theInterpFcn, geogrid_lat, 1, theInterpMethod);
[n, m] = size(grid_x);
FLIPPING = 0;
if FLIPPING
grid_x = flipud(grid_x);
grid_y = flipud(grid_y);
geogrid_lon = flipud(geogrid_lon);
geogrid_lat = flipud(geogrid_lat);
geometry{1} = flipud(geometry{1});
geometry{2} = flipud(geometry{2});
mask = flipud(mask);
imask = flipud(imask);
bathymetry = flipud(bathymetry);
ice_draft = flipud(ice_draft);
ang = flipud(ang);
end
TRANSPOSE = 0;
if TRANSPOSE
grid_x = grid_x';
grid_y = grid_y';
geogrid_lon = geogrid_lon';
geogrid_lat = geogrid_lat';
geometry{1} = geometry{1}';
geometry{2} = geometry{2}';
mask = mask';
imask = imask';
bathymetry = bathymetry';
ice_draft = ice_draft';
ang = ang';
end
xl = max(grid_x(:)) - min(grid_x(:));
el = max(grid_y(:)) - min(grid_y(:));
if ~rem(m, 2), m = m-1; end   % m, n must be odd.
if ~rem(n, 2), n = n-1; end
i_rho = 2:2:m-1; j_rho = 2:2:n-1;
i_psi = 3:2:m-2; j_psi = 3:2:n-2;
i_u   = 3:2:m-2; j_u   = 2:2:n-1;
i_v   = 2:2:m-1; j_v   = 3:2:n-2;
% The xi direction (left-right):
LP = (m-1)/2;   % The rho dimension.
L = LP-1;       % The psi dimension.
% The eta direction (up-down):
MP = (n-1)/2;   % The rho dimension.
M = MP-1;       % The psi dimension.
f = 2.*7.29e-5.*sin(geogrid_lat(j_rho, i_rho).*pi./180);
f(isnan(f))=0;  % doesn't like f to be NaN
% Masking.
%mask = ~~mask;
%land = mask;
%water = ~land;
rmask = mask;
% Calculate other masking arrays.
umask = zeros(size(rmask));
vmask = zeros(size(rmask));
pmask = zeros(size(rmask));
for i = 2:LP
for j = 1:MP
umask(j, i-1) = rmask(j, i) * rmask(j, i-1);
end
end
for i = 1:LP
for j = 2:MP
vmask(j-1, i) = rmask(j, i) * rmask(j-1, i);
end
end
for i = 2:LP
for j = 2:MP
pmask(j-1, i-1) = rmask(j, i) * rmask(j, i-1) * rmask(j-1, i) * rmask(j-1, i-1);
end
end
% Average angle -- We should do this via (x, y) components.
temp = ang;
ang = zeros(n, m);
ang(2:2:end,2:2:end) = temp(2:end,2:end);
gx = geometry{1};   % Spherical distances in meters.
gy = geometry{2};
sx = 0.5*(gx(1:end-1, :) + gx(2:end, :));
sy = 0.5*(gy(:, 1:end-1) + gy(:, 2:end));
pm = 1 ./ sx(:,2:end);
pn = 1 ./ sy(2:end,:);
% pm and pn cannot be Inf, even if on land, so if values
% are Inf, set to an arbitrary non-zero value
pm(isinf(pm))=0.999e-3;
pn(isinf(pn))=0.999e-3;
pm(isnan(pm))=0.999e-3;
pn(isnan(pn))=0.999e-3;
dmde = zeros(size(pm));
dndx = zeros(size(pn));
dmde(2:end-1, :) = 0.5*(1./pm(3:end, :) - 1./pm(1:end-2, :));
dndx(:, 2:end-1) = 0.5*(1./pn(:, 3:end) - 1./pn(:, 1:end-2));
dmde(isinf(dmde))=0;
dndx(isinf(dndx))=0;
dmde(isnan(dmde))=0;
dndx(isnan(dndx))=0;
bb = bathymetry(1:end-1,2:end);
dd = ice_draft(1:end-1,2:end);
% Final check
ii = find((bb+dd).*rmask < min_depth & (bb+dd).*rmask > 0 );
bb(ii) = -dd(ii) + min_depth;
if smooth == 1
    [bb]=smooth_bath(bb,ones(size(bb)),4,rx1in,50);
end
ii = find((bb+dd).*rmask < min_depth & (bb+dd).*rmask > 0 );
bb(ii) = -dd(ii) + min_depth;
%dd(ii) = -bb(ii) + min_depth;

ii = find((bb+dd) < min_depth0);

dd(ii) = -bb(ii) + min_depth0;

ii =  find(dd > 0);
dd(ii) = 0;

bb(ii) = -dd(ii) + min_depth0;





%% SOme check plots:
%if PlotFigs == 1;
%h1 = figure('visible','off'),
%flat_pcolor(x(1:res:end,1:res:end),y(1:res:end,1:res:end),elev_bed(1:res:end,1:res:end));shading flat,
%flat_pcolor(bb)
%colormap(hot_cold_white),
%hold on,
%resg2 = 100;
%plot(Cx(:,1:resg2:end),Cy(:,1:resg2:end),'k');
%hold on, plot(Cx(1:resg2:end,:)',Cy(1:resg2:end,:)','k');
%axis square
%caxis([-1500 1500])
%saveas(h1,'grd_20170809.png','png');
%!cp -p grd_20170809.png ~/public_html/WholeAntarcticModel/.
%end
%% Create netcd file
% btmp = ones(size(bb)).*2000;
% dtmp =  zeros(size(dd));
% dtmp(dd == 0) = 0;
% dtmp(dd < 0) = -0.1;
%
% bb = btmp;
% dd = dtmp;

%% Final write
%htest = ones(size(bathymetry(1:end-1,2:end)')).*900;

GrdName = ['waom' num2str(mr) '_MinDepth' num2str(s.clipping_depths(1)) 'm_rx1' num2str(rx1in) '_grd.nc']
%GrdName = ['waom' num2str(mr) '_MinDepth' num2str(s.clipping_depths(1)) '_grd.nc']

c_grid(LP,MP,GrdName);
nc_write(GrdName,'xl', xl);
nc_write(GrdName,'el', el);
nc_write(GrdName,'f',f');
nc_write(GrdName,'x_rho', grid_x(j_rho, i_rho)');
nc_write(GrdName,'y_rho', grid_y(j_rho, i_rho)');
nc_write(GrdName,'x_psi', grid_x(j_psi, i_psi)');
nc_write(GrdName,'y_psi', grid_y(j_psi, i_psi)');
nc_write(GrdName,'x_u', grid_x(j_u, i_u)');
nc_write(GrdName,'y_u', grid_y(j_u, i_u)');
nc_write(GrdName,'x_v', grid_x(j_v, i_v)');
nc_write(GrdName,'y_v', grid_y(j_v, i_v)');
nc_write(GrdName,'lon_rho', geogrid_lon(j_rho, i_rho)');
nc_write(GrdName,'lat_rho', geogrid_lat(j_rho, i_rho)');
nc_write(GrdName,'lon_psi', geogrid_lon(j_psi, i_psi)');
nc_write(GrdName,'lat_psi', geogrid_lat(j_psi, i_psi)');
nc_write(GrdName,'lon_u', geogrid_lon(j_u, i_u)');
nc_write(GrdName,'lat_u', geogrid_lat(j_u, i_u)');
nc_write(GrdName,'lon_v', geogrid_lon(j_v, i_v)');
nc_write(GrdName,'lat_v', geogrid_lat(j_v, i_v)');
nc_write(GrdName,'h', bb');
%nc_write(GrdName,'h',htest);
nc_write(GrdName,'zice', dd');
nc_write(GrdName,'mask_rho', double(rmask)');
%nc_write(GrdName,'mask_zice', double(imask)');
nc_write(GrdName,'mask_psi', pmask(1:end-1, 1:end-1)');
nc_write(GrdName,'mask_u', umask(1:end, 1:end-1)');
nc_write(GrdName,'mask_v', vmask(1:end-1, 1:end)');
nc_write(GrdName,'angle', ang(j_rho, i_rho)');
nc_write(GrdName,'pm', pm');
nc_write(GrdName,'pn', pn');
nc_write(GrdName,'dmde', dmde');
nc_write(GrdName,'dndx', dndx');

ModelHOME = ['../../../data/preprocessing/interim/']
eval(['!mv ' GrdName ' ' ModelHOME '.']);






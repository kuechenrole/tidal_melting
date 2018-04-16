
addpath ./TMD
addpath ./TMD/FUNCTIONS
addpath ./t_tide_v1.3beta
addpath ./DATA

gfile=[getenv('prodir'),'/waom10_small_grd.nc']
base_date=datenum(2019,1,1);
pred_date=datenum(2019,2,1);
ofile=[getenv('prodir'),'/waom10_tds_2019.nc'];
model_file='DATA/Model_tpxo7.2';
otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,'WAOM10')

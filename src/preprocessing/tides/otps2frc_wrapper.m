addpath(genpath('./tmd_toolbox'))
addpath(genpath('./t_tide'))
proj_dir = getenv('projdir');
prep_dir = [proj_dir,'/data/preprocessing'];

gfile=[prep_dir,'/processed/waom10_grd.nc']
base_date=datenum(2007,01,01);
pred_date=datenum(2007,06,01);
ofile=[prep_dir,'/interim/waom10_tds.nc']
model_file='Model_tpxo7.2'
otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,'WAOM10')

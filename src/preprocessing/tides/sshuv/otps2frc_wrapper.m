addpath(genpath('./tmd_toolbox'))
addpath(genpath('./t_tide'))
proj_dir = getenv('projdir');
prep_dir = [proj_dir,'/data/preprocessing'];

gfile=[prep_dir,'/processed/waom10_small_grd.nc']
base_date=datenum(2007,01,01);
pred_date=datenum(2007,01,15);
ofile=[prep_dir,'/processed/waom10_small_CATS_tds.nc']
%model_file=[getenv('extdir'),'/tpxo/Model_tpxo7.2']
model_file=[getenv('extdir'),'/cats/Model_tpxo']
otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,'WAOM10_small')

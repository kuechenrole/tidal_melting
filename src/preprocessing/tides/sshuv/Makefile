int_dir = $(proj_dir)/data/preprocessing/interim
pro_dir = $(proj_dir)/data/preprocessing/processed

$(pro_dir)/waom5_tds.nc: $(int_dir)/waom5_tds.nc
	ncrename -O -v tide_Ephase,tide_Pphase $^ $@ 
	ncrename -O -v tide_Eamp,tide_Pamp $@

$(int_dir)/waom5_tds.nc: otps2frc_wrapper.m
	matlab -nodesktop -nosplash -r 'run otps2frc_wrapper;quit()'


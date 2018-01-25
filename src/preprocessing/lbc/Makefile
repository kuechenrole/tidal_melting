#pro_dir = $(proj_dir)/data/preprocessing/processed

$(pro_dir)/waom5_bry.nc: make_lbc.m
	matlab -nodesktop -nosplash -r 'run make_lbc; quit()'

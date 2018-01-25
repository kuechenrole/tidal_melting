pro_dir = $(proj_dir)/data/preprocessing/processed

$(pro_dir)/waom5_nudge.nc: make_nudge.py
	python make_nudge.py

$(pro_dir)/waom5_sbc.nc: make_sbc.m
	matlab -nodesktop -nosplash -r 'run make_sbc; quit()'

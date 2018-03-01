proj_dir = getenv('proj_dir');
extern_dir = [proj_dir,'/data/preprocessing/external'];
addpath(genpath('../matlab_tools/'))
Bedmap = '../../../data/preprocessing/external/bedmap2'
addpath(genpath(Bedmap))
RTOPO = '../../../data/preprocessing/external/rtopo/RTopo-2.0.1_30sec_bedrock_topography_S30.nc'

%% Model domain at mesh resolution (mr) in km:
mr =  1
smooth = 1
PlotFigs = 0

make_grid

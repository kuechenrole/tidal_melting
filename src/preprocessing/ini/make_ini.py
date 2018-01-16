import os
import sys
from roms_ini_ecco import roms_ini_ecco
from numpy import array

# set user parameter and call main routine

grd_file = os.path.join(os.environ['proj_dir'],'data','preprocessing','processed','waom10_grd.nc')
out_file = os.path.join(os.environ['proj_dir'],'data','preprocessing','processed','waom10_ini.nc')

print('making ini file for grid file: ' + grd_file)

# Path to ECCO2 files for temperature and salinity 1th January 1995
salt_file = os.path.join(os.environ['proj_dir'],'data','preprocessing','external','ecco2','SALT.nc','SALT.1440x720x50.20070101.nc')
theta_file = os.path.join(os.environ['proj_dir'],'data','preprocessing','external','ecco2','THETA.nc','THETA.1440x720x50.20070101.nc')

# Grid parameters; check grid_file and *.in file to make sure these are correct
Tcline = 50
theta_s = 4.0
theta_b = 0.9
hc = 50
N = 31

# Northernmost index of ECCO2 grid to read (1-based)
nbdry_ecco = 300

# upper and lower bounds for temp and salinity
tempUp = 10
tempLow = -3

saltUp = 34.8
saltLow = 33.2

# define grid bounds of lake vostock
vostock=array([545,325,583,347])

# call function that generates the out file
roms_ini_ecco(grd_file, theta_file, salt_file, out_file, Tcline, theta_s, theta_b, hc, N, nbdry_ecco,(tempLow,tempUp),(saltLow,saltUp),vostock)



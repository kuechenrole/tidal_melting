from netCDF4 import Dataset
from numpy import *
from sys import argv

def add_sponge(grid_path,viscFacMax,diffFacMax,width=5):
    '''Adds a sponge layer frame to an existing grid file with linear increase towards the boundary.
    viscFacMax - maximum factor by which viscosity is multiplied
    diffFacMax - maximum factor by which diffusion is multiplied
    width      - frame width specified in number of cells'''

    # read grid and size
    print('Reading grid.')
    id=Dataset(grid_path,'a')
    lat_rho=id.variables['lat_rho'][:,:]
    num_i=size(lat_rho,0)
    num_j=size(lat_rho,1)



    # set up viscFac and diffFac arrays
    print('Setting up visc and diff Factor arrays.')
    viscFac=ones(shape(lat_rho))
    diffFac=ones(shape(lat_rho))
    for i in range(width):
        for j in range(i,num_j-i):
            viscFac[i,j]=viscFacMax*(1-i/width)
            viscFac[num_i-1-i,j]=viscFacMax*(1-i/width)
            diffFac[i,j]=diffFacMax*(1-i/width)
            diffFac[num_i-1-i,j]=diffFacMax*(1-i/width)
    for j in range(width):
        for i in range(j,num_i-j):
            viscFac[i,j]=viscFacMax*(1-j/width)
            viscFac[i,num_j-1-j]=viscFacMax*(1-j/width)
            diffFac[i,j]=diffFacMax*(1-j/width)
            diffFac[i,num_j-1-j]=diffFacMax*(1-j/width)
                    
    #write visc and diff factor arrays to grid file, update if it already exists 
    print('Writing arrays out to grid file.')
    if 'visc_factor' in id.variables:
        print('Variable "visc_factor" already exists. Updating its value.')
        id.variables['visc_factor'][:,:] = abs(viscFac)
    else:
        id.createVariable('visc_factor', 'f8', ('eta_rho', 'xi_rho'))
        id.variables['visc_factor'][:,:] = abs(viscFac)
    if 'diff_factor' in id.variables:
        print('Variable "diff_factor" already exists. Updating its value.')
        id.variables['diff_factor'][:,:] = abs(diffFac)
    else:
        id.createVariable('diff_factor', 'f8', ('eta_rho', 'xi_rho'))
        id.variables['diff_factor'][:,:] = abs(diffFac)
    id.close()
   

def main():

    script=argv[0]
    grid_path = argv[1]
    viscFacMax = float(argv[2])
    diffFacMax = float(argv[3])
    width = int(argv[4])

    add_sponge(grid_path,viscFacMax,diffFacMax,width)

# commandline interface
if __name__=='__main__':
    main()



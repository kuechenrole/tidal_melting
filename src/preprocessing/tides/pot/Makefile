FC = gfortran
FFLAGS = -O2

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

tidpot: tidpot_load.o tidpot_load_sal.o
	$(FC) $(FFLAGS) -o tidpot tidpot_load.o tidpot_load_sal.o

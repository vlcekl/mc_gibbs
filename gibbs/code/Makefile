FC = gfortran
FFLAGS = -O2 

FILES = rndmod.f90 globmc.f90 molmc.f90 cfgmc.f90 enemc.f90 movemc.f90 measmc.f90
FILSIM = $(FILES) cbmc.f90

#all: cbmc clean

cbmc: $(FILSIM)
	$(FC) $(FFLAGS) -o $@ $(FILSIM)

clean: cbmc
	rm *.mod


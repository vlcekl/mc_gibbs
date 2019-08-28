FC = gfortran
FFLAGS = -O2

FILESLJ = global_mod.f90 model_AM_mod.f90 config_mod.f90 energy_AM_mod.f90 move_AM_mod.f90 measure_mod.f90
FILSILJ = $(FILESLJ) gcmc.f90

all: clean

gcmc_am: $(FILSILJ)
	$(FC) $(FFLAGS) -o $@ $(FILSILJ)

clean: gcmc_am
	rm *.mod


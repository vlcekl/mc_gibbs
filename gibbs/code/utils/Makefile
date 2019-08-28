FC = ifort
FFLAGS = -O3 -warn general

FILES = ../rndmod.f90 ../globmc.f90 ../molmc.f90 ../cfgmc.f90 ../enemc.f90 ../movemc.f90 ../measmc.f90
FILTEST = $(FILES) mcewtest.f90

all: mcewtest clean

mcewtest: $(FILTEST)
	$(FC) $(FFLAGS) -o $@ $(FILTEST)

bf: bf.f90
	$(FC) $@ bf.f90

clean: mcewtest
	rm *.mod


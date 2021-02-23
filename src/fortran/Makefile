# This is an commentary line in a makefile
# Start of the makefile

LIB = ./lib/

# Compiler/Linker settings
FC = gfortran
FCFLAGS = -c -O3 -march=native -fno-underscoring -fPIC #-ffast-math

all: lib

debug: FCFLAGS = -c -Wall -Wextra -Wimplicit-interface -Wno-unused-function -Wno-tabs -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=zero,overflow -finit-real=snan 
#-finit-integer=-9999999
debug: lib

lib: mdle_typedef.o merged-common.o merged-PSTD.o merged-PNLCG.o merged-LBFGS.o merged-PLBFGS.o
	$(FC) -shared -fPIC mdle_typedef.o merged-common.o merged-PSTD.o merged-PNLCG.o merged-LBFGS.o merged-PLBFGS.o -o $(LIB)libOPTIM.so
mdle_typedef.o: mdle_typedef.f90
	$(FC) $(FCFLAGS) mdle_typedef.f90	
merged-common.o: merged-common.f90 
	$(FC) $(FCFLAGS) merged-common.f90
merged-PSTD.o: merged-PSTD.f90 
	$(FC) $(FCFLAGS) merged-PSTD.f90
merged-PNLCG.o: merged-PSTD.o merged-PNLCG.f90 
	$(FC) $(FCFLAGS) merged-PNLCG.f90
merged-LBFGS.o: merged-LBFGS.f90 
	$(FC) $(FCFLAGS) merged-LBFGS.f90
merged-PLBFGS.o: merged-LBFGS.o merged-PLBFGS.f90
	$(FC) $(FCFLAGS) merged-PLBFGS.f90
clean:
	rm  $(LIB)libOPTIM.so *.mod *.o 
# End of the makefile

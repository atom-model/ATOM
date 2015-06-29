# Makefile fuer die Berechnung von PalaeoAtmosphaerenzirkulation
#
# ----------------------------------------------------------------------------


CC = g++
CFLAGS = -g -Wall
LD = g++
#LDFLAGS = -lm -lnetcdf
LDFLAGS = -lm

OBJ = AGCM_main.o Accuracy.o Array.o Array_2D.o Array_1D.o Pressure.o PostProcess_Atmosphere.o Print_Atmosphere.o BC_Atmosphere.o BC_Bathymetry_Atmosphere.o BC_Thermo.o RHS_Atmosphere.o RungeKutta_Atmosphere.o Results_MSL_Atm.o Restore.o MinMax.o



AGCM: $(OBJ)
		$(CC) $(CFLAGS) -o AGCM $(OBJ) $(LDFLAGS)

%.o: %.cpp
		$(CC) $(CFLAGS) -c $<


.PHONY: clean
clean:
		\rm -vrf $(OBJ) AGCM

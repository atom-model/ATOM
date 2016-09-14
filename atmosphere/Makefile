# Makefile fuer die Berechnung von PalaeoAtmosphaerenzirkulation
#
# ----------------------------------------------------------------------------


CC = g++

CFLAGS = -g -Wall
LD = g++
#CFLAGS = -g -Wall -I/usr/local/atom/netcdf-4.4.0/include

# we're just adding all of the potential library locations and hoping that one of them is correct
LDFLAGS = -lm -L/usr/local/lib -L/usr/local/atom/netcdf-4.4.0/lib -lnetcdf

AGCM_OBJ = AGCM_main.o Accuracy_Atm.o Array.o Array_2D.o Array_1D.o Pressure_Atm.o PostProcess_Atm.o Print_Atm.o BC_Atm.o BC_Bath_Atm.o BC_Thermo.o RHS_Atm.o RungeKutta_Atm.o Results_Atm.o Restore_Atm.o MinMax_Atm.o File_NetCDF.o
ATMOSPHERE_OBJ = atmosphere.o Accuracy_Atm.o Array.o Array_2D.o Array_1D.o Pressure_Atm.o PostProcess_Atm.o Print_Atm.o BC_Atm.o BC_Bath_Atm.o BC_Thermo.o RHS_Atm.o RungeKutta_Atm.o Results_Atm.o Restore_Atm.o MinMax_Atm.o File_NetCDF.o AtmosphereConfig.o AtmosphereModel.o

all: AGCM atmosphere

AGCM: $(AGCM_OBJ)
		$(CC) $(CFLAGS) -o AGCM $(AGCM_OBJ) $(LDFLAGS)

atmosphere: $(ATMOSPHERE_OBJ)
		$(CC) $(CFLAGS) -o atmosphere $(ATMOSPHERE_OBJ) $(LDFLAGS)

%.o: %.cpp
		$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
		\rm -vrf $(AGCM_OBJ) $(ATMOSPHERE_OBJ) AGCM atmosphere

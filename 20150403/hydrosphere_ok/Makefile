# Makefile fuer die Berechnung von PalaeoHydrosphaerenzirkulation
#
# ----------------------------------------------------------------------------

CC = g++
CFLAGS = -g -Wall
LD = g++
LDFLAGS = -lm -lnetcdf -lf2c -L/opt/local/lib -I/opt/local/include


OBJ = OGCM_main.o Accuracy.o Array.o Array_2D.o Array_1D.o Pressure.o PostProcess_Hydrosphere.o Print_Hydrosphere.o BC_Hydrosphere.o BC_Bathymetry_Hydrosphere.o IC_Thermohalin.o RHS_Hydrosphere.o RungeKutta_Hydrosphere.o Restore.o MinMax.o Results_MSL_Hyd.o File_NetCDF.o



OGCM: $(OBJ)
		$(CC) $(CFLAGS) -o OGCM $(OBJ) $(LDFLAGS)

%.o: %.cpp
		$(CC) $(CFLAGS) -c $< $(LDFLAGS)


.PHONY: clean
clean:
		\rm -vrf $(OBJ) OGCM

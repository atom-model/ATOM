# Makefile fuer die Berechnung von PalaeoHydrosphaerenzirkulation
#
# ----------------------------------------------------------------------------

CC = g++
CFLAGS = -g -Wall
LD = g++
LDFLAGS = -lm -lnetcdf
#LDFLAGS = -lm


OBJ = OGCM_main.o Accuracy_Hyd.o Array.o Array_2D.o Array_1D.o Pressure_Hyd.o PostProcess_Hyd.o Print_Hyd.o BC_Hyd.o BC_Bath_Hyd.o BC_Thermohalin.o RHS_Hyd.o RungeKutta_Hyd.o Restore_Hyd.o MinMax_Hyd.o Results_Hyd.o File_NetCDF.o

#OBJ = OGCM_main.o Accuracy_Hyd.o Array.o Array_2D.o Array_1D.o Pressure_Hyd.o PostProcess_Hyd.o Print_Hyd.o BC_Hyd.o BC_Bath_Hyd.o BC_Thermohalin.o RHS_Hyd.o RungeKutta_Hyd.o Restore_Hyd.o MinMax_Hyd.o Results_Hyd.o



OGCM: $(OBJ)
		$(CC) $(CFLAGS) -o OGCM $(OBJ) $(LDFLAGS)

%.o: %.cpp
		$(CC) $(CFLAGS) -c $<


.PHONY: clean
clean:
		\rm -vrf $(OBJ) OGCM

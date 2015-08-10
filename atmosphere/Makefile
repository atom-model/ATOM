# Makefile fuer die Berechnung von PalaeoAtmosphaerenzirkulation
#
# ----------------------------------------------------------------------------


CC = g++
CFLAGS = -g -Wall
LD = g++
#LDFLAGS = -lm -lnetcdf
LDFLAGS = -lm

OBJ = AGCM_main.o Accuracy_Atm.o Array.o Array_2D.o Array_1D.o Pressure_Atm.o PostProcess_Atm.o Print_Atm.o BC_Atm.o BC_Bath_Atm.o BC_Thermo.o RHS_Atm.o RungeKutta_Atm.o Results_Atm.o Restore_Atm.o MinMax_Atm.o



AGCM: $(OBJ)
		$(CC) $(CFLAGS) -o AGCM $(OBJ) $(LDFLAGS)

%.o: %.cpp
		$(CC) $(CFLAGS) -c $<


.PHONY: clean
clean:
		\rm -vrf $(OBJ) AGCM

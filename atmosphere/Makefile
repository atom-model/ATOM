# Makefile fuer die Berechnung von PalaeoAtmosphaerenzirkulation
#
# ----------------------------------------------------------------------------


CC = g++

CFLAGS = -g -Wall -I../lib
LD = g++

# we're just adding all of the potential library locations and hoping that one of them is correct
LDFLAGS = -lm -L/usr/local/lib -L/usr/local/atom/netcdf-4.4.0/lib -lnetcdf

COMMON_OBJ = Accuracy_Atm.o ../Array.o Array_2D.o Array_1D.o Pressure_Atm.o PostProcess_Atm.o Print_Atm.o BC_Atm.o BC_Bath_Atm.o BC_Thermo.o RHS_Atm.o RungeKutta_Atm.o Results_Atm.o Restore_Atm.o MinMax_Atm.o File_NetCDF_Atm.o

AGCM_OBJ = AGCM_main.o $(COMMON_OBJ)
ATMOSPHERE_OBJ = atmosphere.o AtmosphereConfig.o AtmosphereModel.o $(COMMON_OBJ)
PYTHON_OBJ = AtmospherePython.o AtmosphereConfig.o AtmosphereModel.o PythonStream.o $(COMMON_OBJ)

all: AGCM atmosphere libagcm.so

AGCM: $(AGCM_OBJ)
		$(CC) $(CFLAGS) -o AGCM $(AGCM_OBJ) $(LDFLAGS)

atmosphere: $(ATMOSPHERE_OBJ)
		$(CC) $(CFLAGS) -o atmosphere $(ATMOSPHERE_OBJ) $(LDFLAGS)

libagcm.so: $(PYTHON_OBJ)
	# below might need ldflags
	$(CC) -shared -Wl -o libagcm.so $(PYTHON_OBJ) $(LDFLAGS)

%.o: %.cpp
		$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
		\rm -vrf $(AGCM_OBJ) $(ATMOSPHERE_OBJ) $(PYTHON_OBJ) AGCM atmosphere agcm.so


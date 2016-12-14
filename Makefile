# Top-level Makefile
# Builds everything: atmosphere, hydrosphere, Python interface and CLI interface

# TODO: only building AGCM right now

CC = g++
CFLAGS = -g -Wall -Ilib -Iatmosphere
LD = g++

# we're just adding all of the potential library locations and hoping that one of them is correct
LDFLAGS = -lm -L/usr/local/lib -L/usr/local/atom/netcdf-4.4.0/lib -lnetcdf

COMMON_OBJ = lib/Array.o lib/Array_2D.o lib/Array_1D.o
ATM_OBJ = atmosphere/Pressure_Atm.o atmosphere/PostProcess_Atm.o atmosphere/Print_Atm.o atmosphere/BC_Atm.o atmosphere/BC_Bath_Atm.o atmosphere/BC_Thermo.o atmosphere/RHS_Atm.o atmosphere/RungeKutta_Atm.o atmosphere/Results_Atm.o atmosphere/Restore_Atm.o atmosphere/MinMax_Atm.o atmosphere/File_NetCDF_Atm.o atmosphere/Accuracy_Atm.o

AGCM_OBJ = atmosphere/AGCM_main.o $(COMMON_OBJ) $(ATM_OBJ)
# ATMOSPHERE_OBJ = atmosphere.o AtmosphereConfig.o AtmosphereModel.o $(COMMON_OBJ)
#PYTHON_OBJ = AtmospherePython.o AtmosphereConfig.o AtmosphereModel.o PythonStream.o $(COMMON_OBJ)

all: AGCM# atmosphere

AGCM: $(AGCM_OBJ)
	$(CC) $(CFLAGS) -o AGCM $(AGCM_OBJ) $(LDFLAGS)

# atm: $(ATMOSPHERE_OBJ)
	# $(CC) $(CFLAGS) -o atm $(ATMOSPHERE_OBJ) $(LDFLAGS)

lib/%.o: lib/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

atmosphere/%.o: atmosphere/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	\rm -vf $(AGCM_OBJ) $(ATMOSPHERE_OBJ) $(PYTHON_OBJ) AGCM atm


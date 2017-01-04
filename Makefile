# Top-level Makefile
# Builds everything: atmosphere, hydrosphere, Python interface and CLI interface

# TODO: don't always enable debugging
CFLAGS = -g -Wall -Ilib -Iatmosphere -Itinyxml2

# we're just adding all of the potential library locations and hoping that one of them is correct
LDFLAGS = -lm -L/usr/local/lib -L/usr/local/atom/netcdf-4.4.0/lib -lnetcdf

# Common files for the shared lib (libatom.a)
LIB_OBJ = lib/Array.o lib/Array_2D.o lib/Array_1D.o
ATM_OBJ = atmosphere/cAtmosphereModel.o atmosphere/Pressure_Atm.o atmosphere/PostProcess_Atm.o atmosphere/Print_Atm.o atmosphere/BC_Atm.o atmosphere/BC_Bath_Atm.o atmosphere/BC_Thermo.o atmosphere/RHS_Atm.o atmosphere/RungeKutta_Atm.o atmosphere/Results_Atm.o atmosphere/Restore_Atm.o atmosphere/MinMax_Atm.o atmosphere/File_NetCDF_Atm.o atmosphere/Accuracy_Atm.o
XML_OBJ = tinyxml2/tinyxml2.o

# Command line objects
CLI_OBJ = cli/atm.o cli/DefaultStream.o
# TODO: add Hydrosphere
# HYD_OBJ = hydrosphere/Accuracy_Hyd.o hydrosphere/BC_Hyd.o File_NetCDF_Hyd.o MinMax_Hyd.o PostProcess_Hyd.o Print_Hyd.o Restore_Hyd.o RungeKutta_Hyd.o BC_Bath_Hyd.o BC_Thermohalin.o IC_Thermohalin.o OGCM_main.o Pressure_Hyd.o RHS_Hyd.o Results_Hyd.o

all: atm python# hyd

libatom.a: $(LIB_OBJ) $(ATM_OBJ) $(XML_OBJ)
	ar rcs libatom.a $(LIB_OBJ) $(ATM_OBJ) $(XML_OBJ)

atm: libatom.a $(CLI_OBJ)
	$(CXX) $(CFLAGS) $(LDFLAGS) -latom -L. -o atm $(CLI_OBJ)

python: libatom.a python/pyatom.so

python/pyatom.so:
	cd python && python setup.py build_ext --inplace








lib/%.o: lib/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

cli/%.o: cli/%.cpp
		$(CXX) $(CFLAGS) -c $< -o $@

atmosphere/%.o: atmosphere/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

tinyxml2/%.o: tinyxml2/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<

.PHONY: clean
clean:
	\rm -vf $(LIB_OBJ) $(ATM_OBJ) $(CLI_OBJ) atm hyd libatom.a
	\rm -vf python/*.so python/*.o python/pyatom.cpp
	\rm -rf python/build/
# AGCM_OBJ = atmosphere/AGCM_main.o $(COMMON_OBJ) $(ATM_OBJ)
# ATMOSPHERE_OBJ = atmosphere.o AtmosphereConfig.o AtmosphereModel.o $(COMMON_OBJ)
#PYTHON_OBJ = AtmospherePython.o AtmosphereConfig.o AtmosphereModel.o PythonStream.o $(COMMON_OBJ)

# all: AGCM# atmosphere

# AGCM: $(AGCM_OBJ)
	# $(CC) $(CFLAGS) -o AGCM $(AGCM_OBJ) $(LDFLAGS)

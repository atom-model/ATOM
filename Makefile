# Top-level Makefile
# Builds everything: atmosphere, hydrosphere, Python interface and CLI interface

# TODO: don't always enable debugging
CFLAGS = -ggdb -Wall -fPIC -std=c++11 -Ilib -Iatmosphere -Ihydrosphere -Itinyxml2

# Common files for the shared lib (libatom.a)
LIB_OBJ = lib/Array.o lib/Array_2D.o lib/Array_1D.o lib/Config.o
ATM_OBJ = atmosphere/cAtmosphereModel.o atmosphere/cAtmosphereConfig.o atmosphere/Pressure_Atm.o atmosphere/PostProcess_Atm.o atmosphere/Print_Atm.o atmosphere/BC_Atm.o atmosphere/BC_Bath_Atm.o atmosphere/BC_Thermo.o atmosphere/RHS_Atm.o atmosphere/RungeKutta_Atm.o atmosphere/Results_Atm.o atmosphere/Restore_Atm.o atmosphere/MinMax_Atm.o atmosphere/Accuracy_Atm.o
HYD_OBJ = hydrosphere/cHydrosphereModel.o hydrosphere/Accuracy_Hyd.o hydrosphere/BC_Hyd.o hydrosphere/MinMax_Hyd.o hydrosphere/PostProcess_Hyd.o hydrosphere/Print_Hyd.o hydrosphere/Restore_Hyd.o hydrosphere/RungeKutta_Hyd.o hydrosphere/BC_Bath_Hyd.o hydrosphere/BC_Thermohalin.o hydrosphere/IC_Thermohalin.o hydrosphere/Pressure_Hyd.o hydrosphere/RHS_Hyd.o hydrosphere/Results_Hyd.o
XML_OBJ = tinyxml2/tinyxml2.o

# Command line objects
ATM_CLI_OBJ = cli/atm.o cli/DefaultStream.o
HYD_CLI_OBJ = cli/hyd.o cli/DefaultStream.o

# Because there are so many parameters, many of the files are autogenerated
PARAM_OUTPUTS = atmosphere/cAtmosphereDefaults.cpp.inc atmosphere/AtmosphereLoadConfig.cpp.inc atmosphere/AtmosphereParams.h.inc hydrosphere/HydrosphereDefaults.cpp.inc hydrosphere/HydrosphereLoadConfig.cpp.inc hydrosphere/HydrosphereParams.h.inc python/atmosphere_params.pxi python/hydrosphere_params.pxi python/atmosphere_pxd.pxi python/hydrosphere_pxd.pxi examples/config.xml python/pyatom.pyx

all: atm hyd python

libatom.a: $(PARAM_OUTPUTS) $(LIB_OBJ) $(ATM_OBJ) $(HYD_OBJ) $(XML_OBJ)
	ar rcs libatom.a $(LIB_OBJ) $(ATM_OBJ) $(HYD_OBJ) $(XML_OBJ)

atm: libatom.a $(ATM_CLI_OBJ)
	$(CXX) $(CFLAGS) $(ATM_CLI_OBJ) -L. -latom $(LDFLAGS) -o atm

# FIXME: this is a bit redundant; there should probably be one binary that does both jobs
hyd: libatom.a $(HYD_CLI_OBJ)
	$(CXX) $(CFLAGS) $(HYD_CLI_OBJ) -L. -latom $(LDFLAGS) -o hyd

$(PARAM_OUTPUTS): param.py
	# explicitly clean dependent files
	rm -f atmosphere/cAtmosphereModel.o atmosphere/cAtmosphereConfig.o hydrosphere/cHydrosphereModel.o hydrosphere/cHydrosphereConfig.o
	python param.py

analyze:
	analyze-build make

python: libatom.a python/pyatom.so

python/pyatom.so: python/pyatom.pyx python/atom.pxd libatom.a
	cd python && python setup.py build_ext --inplace

lib/%.o: lib/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

cli/%.o: cli/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

atmosphere/%.o: atmosphere/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

hydrosphere/%.o: hydrosphere/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

tinyxml2/%.o: tinyxml2/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<

.PHONY: clean
clean:
	\rm -vf $(LIB_OBJ) $(ATM_OBJ) $(HYD_OBJ) $(XML_OBJ) $(ATM_CLI_OBJ) $(HYD_CLI_OBJ) $(PARAM_OUTPUTS) atm hyd libatom.a
	\rm -vf python/*.so python/*.o python/pyatom.cpp
	\rm -rf python/build/

# Given a parameter definition, generates necessary C++, Python and XML bindings
# coding=utf-8


def main():
    # read the input definition

    # name, description, datatype, default
    PARAMS = {
        'common': [
            ('bathymetry_path', '', 'string', 'data/Paleotopography_bathymetry/Golonka_rev210'),
            ('BathymetrySuffix', '', 'string', 'Ma_Golonka.xyz'),
            ('verbose', '', 'bool', False),
            ('output_path', 'directory where model outputs should be placed (must end in /)', 'string', 'output'),
        ],
        'atmosphere': [
            ('velocity_iter_max', '', 'int', 2),
            ('pressure_iter_max', '', 'int', 2),
            ('velocity_iter_max_2D', '', 'int', 2),
            ('pressure_iter_max_2D', '', 'int', 2),
            ('coriolis', 'coriolis force', 'double', 1.0),
            ('centrifugal', 'centrifugal force', 'double', 1.0),
            ('WaterVapour', 'water vapour', 'double', 1.0),
            ('buoyancy', 'buoyancy', 'double', 1.0),
            ('CO2', 'CO2', 'double', 1.0),
            ('epsres', 'accuracy of relative and absolute errors', 'double', 0.00001),
            ('sun', 'while no variable sun position wanted', 'bool', False),
            ('RadiationModel', 'surface temperature computation by a radiation model', 'int', 3),
            ('IceShield', 'compute ice shields? computation of ice shield following the theorie by Milankowitsch', 'bool', False),
            ('declination', 'position of sun axis, today 23,4°, 21.12.: -23,4°, am 21.3. und 23.9.: 0°, 21.6.: +23,4°, in between sin form', 'int', 0),
            ('sun_position_lat', 'position of sun j = 120 means 30°S, j = 60 means 30°N', 'int', 60),
            ('sun_position_lon', 'position of sun k = 180 means 0° or 180° E ( Greenwich, zero meridian )', 'int', 180),
            ('Ma_max', 'parabolic temperature distribution 300 Ma (From Ruddiman)', 'int', 300),
            ('Ma_max_half', 'half of time scale', 'int', 150),
            ('L_atm', 'extension of the atmosphere shell in m, 20000 m / 40 steps = 500 m', 'double', 20000.),
            ('dt', 'time step coincides with the CFL condition', 'double', 0.0001),
            ('dr', 'compares to 500 m hight, 0.0005 * 40 = .02 * 1000 km = 20 km', 'double', 0.0005),
            ('ik', 'solar constant in W/m2', 'double', 1366.),
            ('sigma', 'Stefan-Boltzmann constant W/( m²*K4 )', 'double', 5.670280e-8),
            ('albedo_extra', 'capability of reflection of short wave radiation, global albedo_extra extraterrestric', 'double', 0.15),
            ('epsilon_extra', 'capability of emissions in the atmosphere', 'double', 0.71),
            ('re', 'Reynolds number: ratio viscous to inertia forces, Re = u * L / nue', 'double', 1000.),
            ('ec', 'Eckert number: ratio kinetic energy to enthalpy, Ec = u² / cp T', 'double', 0.00044),
            ('sc_WaterVapour', 'Schmidt number of water vapour, Sc = nue / D', 'double', 0.6),
            ('sc_CO2', 'Schmidt number of CO2', 'double', 0.96),
            ('pr', 'Prandtl number of air for laminar flows', 'double', 0.7179),
            ('g', 'gravitational acceleration of the earth', 'double', 9.8066),
            ('omega', 'rotation number of the earth', 'double', 7.29e-5),
            ('ep', 'ratio of the gas constants of dry air to water vapour [ / ]', 'double', 0.623),
            ('hp', 'water vapour pressure at T = 0°C: E = 6.1 hPa', 'double', 6.1078),
            ('R_Air', 'specific gas constant of air in J/( kg*K ))', 'double', 287.1),
            ('r_air', 'density of dry air in kg/m³ at 20°C', 'double', 1.2041),
            ('R_WaterVapour', 'specific gas constant of water vapour in J/( kg*K ))', 'double', 461.6),
            ('r_water_vapour', 'density of saturated water vapour in kg/m³ at 10°C', 'double', 0.0094),
            ('R_co2', 'specific gas constant of CO2 in J/( kg*4.5K ))', 'double', 188.91),
            ('lv', 'specific latent Evaporation heat ( Condensation heat ) in J/kg', 'double', 2.52e6),
            ('ls', 'specific latent vaporisation heat ( sublimation heat ) in J/kg', 'double', 2.83e6),
            ('cp_l', 'specific heat capacity of dry air at constant pressure and 20°C in J/( kg K )', 'double', 1005.),
            ('Lambda', 'heat transfer coefficient of air in W/m² K )', 'double', 0.0262),
            ('r_water', 'density of water in kg/m³ at 20°C', 'double', 1000.0),
            ('r_co2', 'density of CO2 in kg/m³ at 25°C', 'double', 0.0019767),
            ('gam', 'constant slope of temperature    gam = 0.65 K/100 m', 'double', 0.65),
            ('u_0', 'maximum value of velocity in 15 m/s compares to 54 km/h', 'double', 15.0),
            ('p_0', 'pressure at sea level in hPa', 'double', 1013.25),
            ('t_0', 'temperature in K compare to 0°C', 'double', 273.15),
            ('c_0', 'maximum value of water vapour in kg / kg', 'double', 0.035),
            ('co2_0', 'maximum value of CO2 in ppm at preindustrial times', 'double', 280.0),
            ('ua', 'initial velocity component in r-direction', 'double', 0.0),
            ('va', 'initial velocity component in theta-direction', 'double', 0.0),
            ('wa', 'initial velocity component in phi-direction', 'double', 0.0),
            ('pa', 'initial value for the pressure field', 'double', 0.0),
            ('ca', 'value 1.0 stands for the maximum value of 35 g/kg water vapour', 'double', 0.0),
            ('ta', 'initial value for the temperature field, 1.0 compares to 0° C compares to 273.15 K', 'double', 1.0),
            ('coa', 'initial value of co2 = 1.0 compares to 280 ppm in preindustrial times', 'double', 1.0),
            ('t_cretaceous_max', 'maximum add of mean temperature in °C during cretaceous times', 'double', 10.0),
            ('t_cretaceous', 'value at modern times', 'double', 0.0),
            ('radiation_ocean', 'increase of radiation at equator in W/m²', 'double', 40.0),
            ('radiation_pole', 'negative amount of radiation at poles in W/m²', 'double', -40.0),
            ('radiation_equator', 'positive amount of radiation at equator in W/m²', 'double', 100.0),
            ('t_average', 'mean temperature of the modern earth', 'double', 15.0),
            ('t_equator', 'temperature t_0 = 1.1103 compares to 30.13° C compares to 303.28 K', 'double', 1.1103),
            ('t_pole', 'temperature at the poles t_pole = 0.8 compares to -54.63°C compares to 218.52 K', 'double', 0.8),
            ('t_tropopause', 'temperature in the tropopause, t = 0.78 compares to -60.093°C compares to 213,057 K', 'double', 0.78),
            ('t_land', 'temperature increase on land by 5°C ( 5°C compares to t_land = 0.018305 )', 'double', 0.018305),
            ('c_tropopause', 'minimum water vapour at tropopause c_tropopause = 0.01429 compares to 0.05 g/kg', 'double', 0.0),
            ('c_land', 'water vapour reduction on land ( 90% of the saturation value )', 'double', 0.4),
            ('c_ocean', 'water vapour reduction on sea surface ( 100% of the saturation value )', 'double', 0.5),
            ('co2_average', 'rate of CO2 at preindustrial times', 'double', 372.0),
            ('co2_equator', 'maximum rate of CO2 at sea level at equator, 1. compares to 330 ppm', 'double', 330.0),
            ('co2_tropopause', 'minimum rate CO2 at tropopause 0 ppm', 'double', 0.0),
            ('co2_pole', 'maximum rate of CO2 of the sea surface at poles', 'double', 305.0),
            ('co2_cretaceous', 'value at modern times', 'double', 0.0),
            ('co2_vegetation', 'value compares to 100/600Gt per year on the global surface by vegetation', 'double', 3.0),
            ('co2_ocean', 'value compares to 0.6/600Gt per year on the sea surface', 'double', 0.0),
            ('co2_land', 'value compares to 0.2/600Gt per year on land', 'double', 3.0),
            ('set_sun_position', 'set to true to simulate effect of different sun positions', 'bool', False),
        ],
        'hydrosphere': [
            ('input_path', 'directory where Atmosphere output can be read (must end in /)', 'string', 'output'),
            ('velocity_iter_max', '', 'int', 5),
            ('pressure_iter_max', '', 'int', 5),
            ('velocity_iter_max_2D', '', 'int', 5),
            ('pressure_iter_max_2D', '', 'int', 5),
            ('coriolis', 'computation with Coriolis force', 'double', 1.0),
            ('centrifugal', 'computation with centrifugal force', 'double', 1.0),

        ]
    }

    XML_READ_FUNCS = {
        "string": "FillStringWithElement",
        "double": "FillDoubleWithElement",
        "int": "FillIntWithElement",
        "bool": "FillBoolWithElement"
    }

    def write_cpp_defaults(filename, classname, sections):
        with open(filename, 'w') as f:
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            f.write("void %s::SetDefaultConfig() {\n" % classname)

            for section in sections:
                f.write('\n  // %s section\n' % section)

                for slug, desc, ctype, default in PARAMS[section]:
                    rhs = default
                    if ctype == 'string':
                        rhs = '"%s"' % default
                    elif ctype == 'bool':
                        if default:
                            rhs = 'true'
                        else:
                            rhs = 'false'

                    f.write('  %s = %s;\n' % (slug, rhs))

            f.write("}")

    def write_cpp_load_config(filename, classname, sections):
        with open(filename, 'w') as f:
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")

            for section in sections:
                f.write('\n  // %s section\n' % section)
                element_var_name = 'elem_%s' % section
                f.write('\n  const tinyxml2::XMLElement *%s = atom->FirstChildElement("%s");\n' % (element_var_name, section))
                f.write('\n  if (%s) {\n' % (element_var_name))

                for slug, desc, ctype, default in PARAMS[section]:
                    func_name = XML_READ_FUNCS[ctype]
                    f.write('    Config::%s(%s, "%s", %s);\n' % (func_name, element_var_name, slug, slug))
                f.write("  }\n")

    def write_cpp_headers(filename, sections):
        with open(filename, 'w') as f:
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")

            for section in sections:
                f.write('\n// %s section\n' % section)

                for slug, desc, ctype, default in PARAMS[section]:
                    f.write('%s %s;\n' % (ctype, slug))

    def write_pxi(filename, classname, sections):
        with open(filename, 'w') as f:
            f.write("# THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")

            for section in sections:
                f.write('\n# %s section\n' % section)

                for slug, desc, ctype, default in PARAMS[section]:
                    f.write('property %s\n' % slug)
                    f.write('    def __get__(%s self):\n' % classname)
                    f.write('        self._check_alive()\n')
                    f.write('        return self._thisptr.%s\n' % slug)
                    f.write('\n')
                    f.write('    def __set__(%s self, value):\n' % classname)
                    f.write('        self._check_alive()\n')
                    f.write('        self._thisptr.%s = <%s> value\n' % (slug, ctype))
                    f.write('\n')

    def write_pxd(filename, model, sections):
        with open(filename, 'w') as f:
            # Sadly, Cython docs are incorrect on usage of 'include', so we must include a whole lot of boilerplate
            f.write("""# THIS FILE IS AUTOMATICALLY GENERATED BY param.py
# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME
cdef extern from "c%sModel.h":
    cppclass c%sModel:
        c%sModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(const char *filename)
        void Run()
        void RunTimeSlice(int time_slice)
""" % (model, model, model))

            for section in sections:
                f.write('        # %s section\n' % section)

                for slug, desc, ctype, default in PARAMS[section]:
                    f.write('        %s %s\n' % (ctype, slug))

    def write_config_xml(filename):
        with open(filename, 'w') as f:
            f.write("""<!-- THIS FILE IS GENERATED AUTOMATICALLY BY param.py. DO NOT EDIT. -->
<atom>
""")

            for section in sections:
                f.write('    <%s>\n' % section)

                for slug, desc, ctype, default in PARAMS[section]:
                    if ctype == 'bool':
                        default = str(default).lower()  # Python uses True/False, C++, uses true/false

                    f.write('        <%s>%s</%s>  <!-- %s (%s) -->\n' % (slug, default, slug, desc, ctype))

                f.write('    </%s>\n' % section)
            f.write('</atom>')

    atmosphere_sections = ['common', 'atmosphere']
    hydrosphere_sections = ['common', 'hydrosphere']

    for filename, classname, sections in [
        ('atmosphere/cAtmosphereDefaults.cpp.inc', 'cAtmosphereModel', atmosphere_sections),
        ('hydrosphere/cHydrosphereDefaults.cpp.inc', 'cHydrosphereModel', hydrosphere_sections)
    ]:
        write_cpp_defaults(filename, classname, sections)

    for filename, classname, sections in [
        ('atmosphere/AtmosphereLoadConfig.cpp.inc', 'cAtmosphereModel', atmosphere_sections),
        ('hydrosphere/HydrosphereLoadConfig.cpp.inc', 'cHydrosphereModel', hydrosphere_sections)
    ]:
        write_cpp_load_config(filename, classname, sections)

    for filename, sections in [
        ('atmosphere/AtmosphereParams.h.inc', atmosphere_sections),
        ('hydrosphere/HydrosphereParams.h.inc', hydrosphere_sections)
    ]:
        write_cpp_headers(filename, sections)

    for filename, classname, sections in [
        ('python/atmosphere_params.pxi', 'Atmosphere', atmosphere_sections),
        ('python/hydrosphere_params.pxi', 'Hydrosphere', hydrosphere_sections)
    ]:
        write_pxi(filename, classname, sections)

    for filename, model, sections in [
        ('python/atmosphere_pxd.pxi', 'Atmosphere', atmosphere_sections),
        ('python/hydrosphere_pxd.pxi', 'Hydrosphere', hydrosphere_sections)
    ]:
        write_pxd(filename, model, sections)

    write_config_xml('examples/config.xml')

if __name__ == '__main__':
    main()

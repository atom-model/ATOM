from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension("pyatom",
              [
                  'pyatom.pyx',
                  'cAtmosphereModel.cpp',
                  '../lib/Array.cpp',
                  '../lib/Array_2D.cpp',
                  '../lib/Array_1D.cpp',
                  '../atmosphere/Accuracy_Atm.cpp',
                  '../atmosphere/File_NetCDF_Atm.cpp',
                  '../atmosphere/Pressure_Atm.cpp',
                  '../atmosphere/PostProcess_Atm.cpp',
                  '../atmosphere/Print_Atm.cpp',
                  '../atmosphere/BC_Atm.cpp',
                  '../atmosphere/BC_Bath_Atm.cpp',
                  '../atmosphere/BC_Thermo.cpp',
                  '../atmosphere/RHS_Atm.cpp',
                  '../atmosphere/RungeKutta_Atm.cpp',
                  '../atmosphere/Results_Atm.cpp',
                  '../atmosphere/Restore_Atm.cpp',
                  '../atmosphere/MinMax_Atm.cpp',
              ],
              language='c++',
              libraries=['netcdf'],
              include_dirs=['../atmosphere', '../lib']
              )]

setup(
    name='pyatom',
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)

# TODO: for distribution, include pyatom.cpp
# http://cython.readthedocs.io/en/latest/src/reference/compilation.html#distributing-cython-modules
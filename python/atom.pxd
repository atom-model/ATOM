from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "cAtmosphereModel.h":
    cppclass cAtmosphereModel:
        cAtmosphereModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(const char *filename)
        void Run()
        void RunTimeSlice(int time_slice)

include "atmosphere_pxd.pxi"


include "hydrosphere_pxd.pxi"

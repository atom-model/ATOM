from libcpp.string cimport string

cdef extern from "cAtmosphereModel.h":
    cppclass cAtmosphereModel:
        cAtmosphereModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(string& filename)
        void Run()

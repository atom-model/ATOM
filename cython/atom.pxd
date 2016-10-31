# Using a .pxd file gives us a separate namespace for
# the C++ declarations. Using a .pxd file also allows
# us to reuse the declaration in multiple .pyx modules.

from libcpp.string cimport string

cdef extern from "cAtmosphereModel.h":
    cppclass cAtmosphereModel:
        # int x,y
        cAtmosphereModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(string& filename)
        void Run()
        # int Multiply(int a, int b)

# distutils: language = c++
# distutils: sources = cAtmosphereModel.cpp
# distutils: include_dirs = ../atmosphere ../lib

import sys

from atom cimport cAtmosphereModel, cHydrosphereModel

from libcpp.string cimport string

class Model:
    """
    ATOM Model Object
    """
    def __init__(self):
        self.atmosphere = Atmosphere()
        self.hydrosphere = Hydrosphere()

    def run(self, t, output_path):
        self.atmosphere.output_path = output_path
        self.atmosphere.run_time_slice(t)

        self.hydrosphere.input_path = output_path
        self.hydrosphere.output_path = output_path
        self.hydrosphere.bathymetry_path = '../data/Paleotopography_bathymetry/Golonka_rev210/'  #  FIXME URGENT read this from XML config
        self.hydrosphere.run_time_slice(t)

    def load_config(self, file_name):
        self.atmosphere.load_config(file_name)
        self.hydrosphere.load_config(file_name)

cdef class Atmosphere:
    """ 
    Cython wrapper class for C++ class Atmosphere
    """
    cdef:
        cAtmosphereModel *_thisptr

    def __cinit__(Atmosphere self):
        # Initialize the "this pointer" to NULL so __dealloc__
        # knows if there is something to deallocate. Do not 
        # call new TestClass() here.
        self._thisptr = NULL
        
    def __init__(Atmosphere self):
        # Constructing the C++ object might raise std::bad_alloc
        # which is automatically converted to a Python MemoryError
        # by Cython. We therefore need to call "new TestClass()" in
        # __init__ instead of __cinit__.
        self._thisptr = new cAtmosphereModel() 

    def __dealloc__(Atmosphere self):
        # Only call del if the C++ object is alive, 
        # or we will get a segfault.
        if self._thisptr != NULL:
            del self._thisptr
            
    cdef int _check_alive(Atmosphere self) except -1:
        # Beacuse of the context manager protocol, the C++ object
        # might die before PyAtmosphere self is reclaimed.
        # We therefore need a small utility to check for the
        # availability of self._thisptr
        if self._thisptr == NULL:
            raise RuntimeError("Wrapped C++ object is deleted")
        else:
            return 0

    def load_config(Atmosphere self, str filename):
        self._check_alive()
        self._thisptr.LoadConfig(filename)
        return None
    
    def run(Atmosphere self):
        self._check_alive()
        self._thisptr.Run()
        return None

    def run_time_slice(Atmosphere self, int t):
        self._check_alive()
        self._thisptr.RunTimeSlice(t)
        return None

    # include "atmosphere_params.pxi"

    # The context manager protocol allows us to precisely
    # control the liftetime of the wrapped C++ object. del
    # is called deterministically and independently of 
    # the Python garbage collection.

    def __enter__(Atmosphere self):
        self._check_alive()
        return self
    
    def __exit__(Atmosphere self, exc_tp, exc_val, exc_tb):
        if self._thisptr != NULL:
            del self._thisptr 
            self._thisptr = NULL # inform __dealloc__
        return False # propagate exceptions

cdef class Hydrosphere:
    """ 
    Cython wrapper class for C++ class Hydrosphere
    """
    cdef:
        cHydrosphereModel *_thisptr

    def __cinit__(Hydrosphere self):
        # Initialize the "this pointer" to NULL so __dealloc__
        # knows if there is something to deallocate. Do not 
        # call new TestClass() here.
        self._thisptr = NULL
        
    def __init__(Hydrosphere self):
        # Constructing the C++ object might raise std::bad_alloc
        # which is automatically converted to a Python MemoryError
        # by Cython. We therefore need to call "new TestClass()" in
        # __init__ instead of __cinit__.
        self._thisptr = new cHydrosphereModel() 

    def __dealloc__(Hydrosphere self):
        # Only call del if the C++ object is alive, 
        # or we will get a segfault.
        if self._thisptr != NULL:
            del self._thisptr
            
    cdef int _check_alive(Hydrosphere self) except -1:
        # Beacuse of the context manager protocol, the C++ object
        # might die before PyHydrosphere self is reclaimed.
        # We therefore need a small utility to check for the
        # availability of self._thisptr
        if self._thisptr == NULL:
            raise RuntimeError("Wrapped C++ object is deleted")
        else:
            return 0

    def load_config(Hydrosphere self, str filename):
        self._check_alive()
        self._thisptr.LoadConfig(filename)
        return None
    
    def run(Hydrosphere self):
        self._check_alive()
        self._thisptr.Run()
        return None

    def run_time_slice(Hydrosphere self, int t):
        self._check_alive()
        self._thisptr.RunTimeSlice(t)
        return None

    # include "hydrosphere_params.pxi"

    # The context manager protocol allows us to precisely
    # control the liftetime of the wrapped C++ object. del
    # is called deterministically and independently of 
    # the Python garbage collection.

    def __enter__(Hydrosphere self):
        self._check_alive()
        return self
    
    def __exit__(Hydrosphere self, exc_tp, exc_val, exc_tb):
        if self._thisptr != NULL:
            del self._thisptr 
            self._thisptr = NULL # inform __dealloc__
        return False # propagate exceptions


cdef public void PythonPrint(const char *s):
    print s,
    sys.stdout.flush()

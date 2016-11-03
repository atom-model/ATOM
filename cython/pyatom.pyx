# distutils: language = c++
# distutils: sources = cAtmosphereModel.cpp
# distutils: include_dirs = ../atmosphere ../lib

import sys

from atom cimport cAtmosphereModel

from libcpp.string cimport string


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

    '''
    property x:
    
        # Here we use a property to expose the public member
        # x of TestClass to Python
        
        def __get__(Atmosphere self):
            self._check_alive()
            return self._thisptr.x
        
        def __set__(PyAtmosphere self, value):
            self._check_alive()
            self._thisptr.x = <int> value

    property y:
    
        # Here we use a property to expose the public member
        # y of TestClass to Python
    
        def __get__(PyAtmosphere self):
            self._check_alive()
            return self._thisptr.y
    
        def __set__(PyAtmosphere self, value):
            self._check_alive()
            self._thisptr.y = <int> value

            
    def Multiply(Atmosphere self, int a, int b):
        self._check_alive()
        return self._thisptr.Multiply(a,b)
            
    '''
    def load_config(Atmosphere self, str filename):
        self._check_alive()
        # cdef string cfilename = filename

        print(filename)

        self._thisptr.LoadConfig(filename)

        return None
    
    def run(Atmosphere self):
        self._check_alive()
        self._thisptr.Run()
        return None

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

cdef public void PythonPrint(const char *s):
    print s,
    sys.stdout.flush()

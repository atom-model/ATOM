// Python version: all output goes through Python

#include <iostream>
#include <fstream>

#include <Python.h>
#include "PythonStream.h"
#include "pyatom.h"

using namespace std;

void PythonStream::OverrideCout() {
    // FIXME: leaks a PythonStream object
    // also permanently mangles the output, but that's fine in this application
    cout.rdbuf(new PythonStream());
}

int PythonStream::sync() {
    // call the python print function with this->str()
    PythonPrint(this->str().c_str());

    return 0;
}

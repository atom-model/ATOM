// cli version: does nothing (cout goes to cout)

#include "PythonStream.h"

void PythonStream::OverrideCout() { }

int PythonStream::sync() {
    return 0;
}

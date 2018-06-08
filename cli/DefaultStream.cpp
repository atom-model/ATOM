// cli version: does nothing (cout goes to cout)

#include "PythonStream.h"

bool PythonStream::is_enable(){
    return false;
}

void PythonStream::OverrideCout() { }

int PythonStream::sync() {
    return 0;
}

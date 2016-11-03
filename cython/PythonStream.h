// This header is included in all compilations, but the Python/cli version gets
// a different cpp implementation.

#ifndef PYTHONSTREAM_H
#define PYTHONSTREAM_H

#include <sstream>

class PythonStream : public std::stringbuf
{
public:
    static void OverrideCout();

private:
    virtual int sync();
};

#endif

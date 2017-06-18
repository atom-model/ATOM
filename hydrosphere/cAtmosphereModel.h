#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include <string>

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class cAtmosphereModel {
public:
    const char *filename;

    cAtmosphereModel();
    ~cAtmosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    #include "AtmosphereParams.h.inc"

private:
    void SetDefaultConfig();
};

#endif

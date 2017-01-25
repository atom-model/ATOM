#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include <string>

#include "../tinyxml2/tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class cAtmosphereModel {
public:
    cAtmosphereModel();
    ~cAtmosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    #include "AtmosphereParams.h.inc"

    // TODO: j_sun - priority summer vs winter parameter
    string output_path;

private:
    void SetDefaultConfig();
};

#endif

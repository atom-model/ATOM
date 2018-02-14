#ifndef CHYDROSPHEREMODEL_H
#define CHYDROSPHEREMODEL_H

#include <string>

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class cHydrosphereModel {
public:
    const char *filename;

    cHydrosphereModel();
    ~cHydrosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    #include "HydrosphereParams.h.inc"

private:
    void SetDefaultConfig();
};

#endif

#ifndef CHYDROSPHEREMODEL_H
#define CHYDROSPHEREMODEL_H

#include "../tinyxml2/tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class cHydrosphereModel {
public:
    cHydrosphereModel();
    ~cHydrosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int Ma);

    // CONFIGURATION
    string output_path;
    string input_path;
    string bathymetry_path;
    //TODO

    // SIMULATION PARAMETERS
    // TODO

    // PHYSICAL PARAMETERS
    // To modify the defaults, see cHydrosphereModel.cpp
    // TODO
};

#endif

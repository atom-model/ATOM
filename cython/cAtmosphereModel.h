#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

using namespace std;

class cAtmosphereModel {
public:
    cAtmosphereModel();
    ~cAtmosphereModel();
    
    // FUNCTIONS
    void LoadConfig(const string& filename);
    void Run();

    // CONFIGURATION
    string inputPath;
    string outputPath;

    // parameters
    // To modify the defaults, see cAtmosphereModel.cpp
    double coriolis;                           // computation with Coriolis force
    double centrifugal;                        // computation with centrifugal force
    double WaterVapour;                        // computation with water vapour
    double buoyancy;                           // computation with buoyancy
    double CO2;                                // computation with CO2

    // TODO: j_sun - priority summer vs winter parameter
};

#endif

#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

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
    bool verbose;

    // SIMULATION PARAMETERS
    int velocity_iter_max;
    int pressure_iter_max;

    // PHYSICAL PARAMETERS
    // To modify the defaults, see cAtmosphereModel.cpp
    double coriolis;
    double centrifugal;
    double WaterVapour;
    double buoyancy;
    double CO2;

    // TODO: j_sun - priority summer vs winter parameter

private:
    void FillBoolWithElement(const XMLElement *parent, const char *name, bool &dest) const;
    void FillDoubleWithElement(const XMLElement *parent, const char *name, double &dest) const;
    void FillIntWithElement(const XMLElement *parent, const char *name, int &dest) const;
};

#endif

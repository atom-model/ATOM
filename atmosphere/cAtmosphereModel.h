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
    void WriteConfig(const char *filename) const;
    void Run();
    void RunTimeSlice(int time_slice);

    // CONFIGURATION
    bool verbose;

    string bathymetry_path;
    string bathymetry_suffix;
    string modern_bathymetry_file;

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
    string output_path;

private:
    void SetDefaultConfig();
    
    void FillBoolWithElement(const XMLElement *parent, const char *name, bool &dest) const;
    void FillDoubleWithElement(const XMLElement *parent, const char *name, double &dest) const;
    void FillIntWithElement(const XMLElement *parent, const char *name, int &dest) const;
    void FillStringWithElement(const XMLElement *parent, const char *name, string &dest) const;

};

#endif

#ifndef ATMOSPHERE_CONFIG_H
#define ATMOSPHERE_CONFIG_H

using namespace std;

class AtmosphereConfig
{
    public:
        AtmosphereConfig();
        void LoadConfigFromFile(const string &);

        // parameters
        // To modify the defaults, see AtmosphereConfig.cpp
        double coriolis;                                   // computation with Coriolis force
        double centrifugal;                            // computation with centrifugal force
        double WaterVapour;                        // computation with water vapour
        double buoyancy;                               // computation with buoyancy
        double CO2;                                        // computation with CO2
};

#endif

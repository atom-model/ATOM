#include <iostream>
#include <string>

#include "AtmosphereConfig.h"

using namespace std;

AtmosphereConfig::AtmosphereConfig() {
    coriolis = 1.;                                   // computation with Coriolis force
    centrifugal = 1.;                            // computation with centrifugal force
    WaterVapour = 1.;                        // computation with water vapour
    buoyancy = 1.;                               // computation with buoyancy
    CO2 = 1.;                                        // computation with CO2
}

void AtmosphereConfig::LoadConfigFromFile(const string& filename) {
    cout << "TODO STUB loadconfig " << filename << endl;
}

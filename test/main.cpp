#include <iostream>
#include "cAtmosphereModel.h"
#include "BC_Thermo.h"

class AtomTest
{
public:
    void run()
    {
        cout << "Hello test" << endl;
        cAtmosphereModel model;
        model.LoadConfig("../benchmark/config_atm.xml");

        BC_Thermo thermo(model);
        std::cout.precision(10);
        std::cout <<"the increment: "<< thermo.get_temperature_increment(10) << std::endl;
    }
};

int main(int argc, char **argv) {
    AtomTest test;
    test.run();
}


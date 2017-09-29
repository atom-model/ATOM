#include <iostream>
#include <stdlib.h>

#include "cHydrosphereModel.h"

int main(int argc, char **argv) {
    cHydrosphereModel model;

    if ( argc != 2 )
    {
        std::cout << std::endl << "ATOM Hydrosphere Model" << std::endl;
        std::cout << std::endl;
        std::cout << "Invalid Command Line Parameter" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "\t" << "./hyd <<XML configuration file path>>" << std::endl;
        std::cout << "\t" << "For example: ./hyd config_hyd.xml" << std::endl;
        std::cout << std::endl;
        exit ( 1 );
    }

    model.LoadConfig(argv[1]);
    model.Run();
}


#include <iostream>
#include <stdlib.h>

#include "cAtmosphereModel.h"

int main(int argc, char **argv) 
{
    cAtmosphereModel model;

    if ( argc != 2 )
    {
        std::cout << std::endl << "ATOM Atmosphere Model" << std::endl;
        std::cout << std::endl;
        std::cout << "Invalid Command Line Parameter" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "\t" << "./atm <<XML configuration file path>>" << std::endl;
        std::cout << "\t" << "For example: ./atm config_atm.xml" << std::endl;
        std::cout << std::endl;
        exit ( 1 );
	}
    try{    
        model.LoadConfig(argv[1]);
        model.Run();
    }catch(const std::exception &exc){
        std::cerr << exc.what() << std::endl;
    }
}

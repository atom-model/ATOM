#include <iostream>
#include <stdlib.h>

#include "cHydrosphereModel.h"

int main(int argc, char **argv) {
    cHydrosphereModel model;

    if ( argc < 2 )
    {
        std::cout << endl;
        std::cout << "missing command line parameters\n";
        std::cout << endl;
        exit ( 1 );
	}

    std::cout << endl;
    std::cout << "ATOM hydrosphere model\n";
    std::cout << "\n";
    std::cout << "Usage:\n";
    std::cout << "\t" << argv[0] << "                      " << " <program name>\n";
    std::cout << "\t" << argv[1] << "                  " << " <XML configuration path>\n";
    std::cout << "\t" << argv[2] << "    " << " <XML file name>\n" << "\n";
    std::cout << endl;

    model.LoadConfig(argv[2]);
    model.Run();
}

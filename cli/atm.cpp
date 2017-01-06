#include <iostream>

#include "cAtmosphereModel.h"

int main(int argc, char **argv) {
    cAtmosphereModel model;

    if (argc != 2) {
        std::cout << "ATOM atmosphere model\n";
        std::cout << "\n";
        std::cout << "Usage:\n";
        std::cout << "\t" << argv[0] << " <XML configuration path>\n";
        std::cout << "\n";
        return 1;
    }

    std::cout << "Loading configuration from " << argv[1] << "\n";

    model.LoadConfig(argv[1]);
    model.Run();
}

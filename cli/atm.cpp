#include <iostream>

#include "cAtmosphereModel.h"

/*
void usage() {
    std::cout << "with no command line args, behaves exactly like the old one\n";
    std::cout << "can also take a config file as an argument\n";
}
*/

int main(int argc, char **argv) {
    cAtmosphereModel model;

    std::cout << "ATOM atmosphere model\n";
    std::cout << "Usage:\n";
    std::cout << "\t%s [XML configuration path]\n";
    std::cout << "If XML path is not specified, default parameters will be used\n";

    if (argc == 2) {
        model.LoadConfig(argv[1]);
    }

    model.Run();
}

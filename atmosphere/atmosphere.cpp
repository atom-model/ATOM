#include <iostream>
#include "AtmosphereConfig.h"
#include "AtmosphereModel.h"

using namespace std;

void usage() {
    cout << "with no command line args, behaves exactly like the old one" << endl;
    cout << "can also take a config file as an argument" << endl;
}

int main(int argc, char **argv) {
    AtmosphereConfig config;

    if (argc == 2) {
        config.LoadConfigFromFile(argv[1]);
    }

    AtmosphereModel(config);
}

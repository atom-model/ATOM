#include <iostream>
#include <stdexcept>

#include "cAtmosphereModel.h"
#include "Config.h"
#include "tinyxml2.h"

void cAtmosphereModel::LoadConfig(const char *filename) {
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    cout << "err " << err << "\n";
    if (err) {
        doc.PrintError();
        throw std::invalid_argument("couldn't load config file");
    }

    XMLElement *atom = doc.FirstChildElement("atom");
    if (!atom) {
        return;
    }

    #include "AtmosphereLoadConfig.cpp.inc"
}

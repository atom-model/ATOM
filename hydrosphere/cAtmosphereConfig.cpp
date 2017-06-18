#include <iostream>
#include <stdexcept>

#include "cAtmosphereModel.h"
#include "Config.h"
#include "tinyxml2.h"

void cAtmosphereModel::LoadConfig(const char *filename) {
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    cout << "reading config-xml-file in cAtmosphereModel::LoadConfig ======= error message = err " << err << "\n\n";
    if (err) {
        doc.PrintError();
        throw std::invalid_argument("couldn't load config file");
    }

    XMLElement *atom = doc.FirstChildElement("atom");
    if (!atom) {
        return;
    }

	XMLElement* elem_common = doc.FirstChildElement( "atom" )->FirstChildElement( "elem_common" );
	if (!elem_common) {
		return;
	}

	XMLElement* elem_atmosphere = doc.FirstChildElement( "elem_common" )->FirstChildElement( "elem_atmosphere" );
	if (!elem_atmosphere) {
		return;
	}


	#include "AtmosphereLoadConfig.cpp.inc"
}

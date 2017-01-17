#include <iostream>
#include <stdexcept>

#include "cAtmosphereModel.h"
#include "tinyxml2.h"

void cAtmosphereModel::SetDefaultConfig() {
    // common parameters
    verbose = false;

    bathymetry_path = "data/Paleotopography_bathymetry/Golonka_rev210";
    bathymetry_suffix = "Ma_Golonka.xyz";

    // simulation parameters
    velocity_iter_max = 2;
    pressure_iter_max = 2;

    // physical parameters
    coriolis = 1.;
    centrifugal = 1.;
    WaterVapour = 1.;
    buoyancy = 1.;
    CO2 = 1.;
}

void cAtmosphereModel::WriteConfig(const char *filename) const {
    std::cout << "Writing config to " << filename << "\n";

    XMLDocument doc;

    XMLElement *atom_element = doc.NewElement("atom");
    doc.InsertEndChild(atom_element);

    // TODO common parameters

    XMLElement *atmosphere_element = doc.NewElement("atmosphere");
    atom_element->InsertEndChild(atmosphere_element);

    // TODO: automatically generate the list of parameters and types
    XMLElement *elem = doc.NewElement("coriolis");
    atmosphere_element->InsertEndChild(elem);
    elem->SetText(coriolis);

    XMLError err = doc.SaveFile(filename);
    if (err) {
        // TODO: proper error handling
        cout << "error " << err << " while writing config\n";
        abort();
    }
}

void cAtmosphereModel::FillDoubleWithElement(const XMLElement *parent, const char *name, double &dest) const {
    const XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (text) {
        dest = stod(text);
    }
}

void cAtmosphereModel::FillIntWithElement(const XMLElement *parent, const char *name, int &dest) const {
    const XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (text) {
        dest = stoi(text);
    }
}

void cAtmosphereModel::FillStringWithElement(const XMLElement *parent, const char *name, string &dest) const {
    const XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (text) {
        dest = text;
    }
}

void cAtmosphereModel::FillBoolWithElement(const XMLElement *parent, const char *name, bool &dest) const {
    const XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (0 == strcmp(text, "false")) {
        dest = false;
    } else if (0 == strcmp(text, "true")) {
        dest = true;
    } else {
        cout << "ERROR: unknown value '" << text << "' for config item '" << name << "'; I expect 'true' or 'false'" << endl;
    }
}

void cAtmosphereModel::LoadConfig(const char *filename) {
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    cout << "err " << err << "\n";
    if (err) {
        doc.PrintError();
        throw std::invalid_argument("couldn't load config file");
    }

    // TODO you should catch exceptions:
    // libc++abi.dylib: terminating with uncaught exception of type std::invalid_argument: stod: no conversion

    XMLElement *atom = doc.FirstChildElement("atom");

    // SIMULATION PARAMETERS
    XMLElement *common = atom->FirstChildElement("common");

    FillStringWithElement(common, "BathymetryPath", bathymetry_path);
    FillStringWithElement(common, "BathymetrySuffix", bathymetry_suffix);
    FillStringWithElement(common, "ModernBathymetryFile", modern_bathymetry_file);

    // ATMOSPHERE PHYSICAL PARAMETERS
    XMLElement *atm = atom->FirstChildElement("atmosphere");

    FillIntWithElement(atm, "velocity_iter_max", velocity_iter_max);
    FillIntWithElement(atm, "pressure_iter_max", pressure_iter_max);

    FillDoubleWithElement(atm, "coriolis", coriolis);
    FillDoubleWithElement(atm, "centrifugal", centrifugal);
    FillDoubleWithElement(atm, "WaterVapour", WaterVapour);
    FillDoubleWithElement(atm, "buoyancy", buoyancy);
    FillDoubleWithElement(atm, "CO2", CO2);
}


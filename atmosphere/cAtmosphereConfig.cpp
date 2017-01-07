#include <iostream>

#include "cAtmosphereModel.h"
#include "tinyxml2.h"

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

#include "Config.h"

#include <iostream>
#include <stdexcept>

void Config::FillDoubleWithElement(const tinyxml2::XMLElement *parent, const char *name, double &dest) {
    const tinyxml2::XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (text) {
        try {
            dest = std::stod(text);
        } catch (std::invalid_argument) {
            std::cout << "ERROR: while reading XML file, could not convert '" << text << "' to double for parameter " << name << "\n";
            std::exit(1);
        }
    }
}

void Config::FillIntWithElement(const tinyxml2::XMLElement *parent, const char *name, int &dest) {
    const tinyxml2::XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (text) {
        try {
            dest = std::stoi(text);
        } catch (std::invalid_argument) {
            std::cout << "ERROR: while reading XML file, could not convert '" << text << "' to int for parameter " << name << "\n";
            std::exit(1);
        }
    }
}

void Config::FillStringWithElement(const tinyxml2::XMLElement *parent, const char *name, std::string &dest) {
    const tinyxml2::XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (text) {
        dest = text;
    }
}

void Config::FillBoolWithElement(const tinyxml2::XMLElement *parent, const char *name, bool &dest) {
    const tinyxml2::XMLElement *elem = parent->FirstChildElement(name);
    if (!elem) {
        return;
    }

    const char *text = elem->GetText();
    if (0 == strcmp(text, "false")) {
        dest = false;
    } else if (0 == strcmp(text, "true")) {
        dest = true;
    } else {
        std::cout << "ERROR: while reading XML file, could not convert '" << text << "' to bool (true/false) for parameter " << name << "\n";
        std::exit(1);
    }
}

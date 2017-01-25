#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>

#include "tinyxml2.h"

class Config {
public:
    static void FillBoolWithElement(const tinyxml2::XMLElement *parent, const char *name, bool &dest);
    static void FillDoubleWithElement(const tinyxml2::XMLElement *parent, const char *name, double &dest);
    static void FillIntWithElement(const tinyxml2::XMLElement *parent, const char *name, int &dest);
    static void FillStringWithElement(const tinyxml2::XMLElement *parent, const char *name, std::string &dest);
};

#endif /* CONFIG_H_ */

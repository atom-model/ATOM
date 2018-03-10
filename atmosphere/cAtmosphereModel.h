#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include <string>
#include <vector>

#include "tinyxml2.h"
#include "Array.h"

using namespace std;
using namespace tinyxml2;

class cAtmosphereModel {
public:
    const char *filename;

    cAtmosphereModel();
    ~cAtmosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    #include "AtmosphereParams.h.inc"

private:
    void SetDefaultConfig();

    std::vector<std::vector<double> > m_node_weights;
    float GetMeanTemperature(int jm, int km, Array &t);
    void CalculateNodeWeights(int jm, int km);
};

#endif

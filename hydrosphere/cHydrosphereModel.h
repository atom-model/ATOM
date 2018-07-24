#ifndef CHYDROSPHEREMODEL_H
#define CHYDROSPHEREMODEL_H

#include <string>
#include <vector>

#include "Array.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class cHydrosphereModel {
public:

    cHydrosphereModel();
    ~cHydrosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    #include "HydrosphereParams.h.inc"

private:
    void SetDefaultConfig();

    //code for future, do not delete
    //std::vector<Array*> old_arrays_3d, new_arrays_3d, old_arrays_2d, new_arrays_2d;
};

#endif

#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include <string>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "tinyxml2.h"

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
    void PrintMaxMinValues();
    #include "AtmosphereParams.h.inc"

private:
    void SetDefaultConfig();

    int im, jm, km, nm;
    int j_res, k_res;    

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr, dt;
    static const double the0, phi0, r0;
    double coeff_mmWS;
    int *im_tropopause;  // location of the tropopause

    Array_2D albedo; // albedo = reflectivity
};

#endif

/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to surveil the accuracy of the iterations
*/
#include <tuple>

#include "Array.h"
#include "Array_1D.h"

#ifndef _ACCURACY_
#define _ACCURACY_

using namespace std;
class Accuracy_Atm{
    private:
        int im, jm, km;
        int i_res, j_res, k_res;
        double dr, dthe, dphi;
        double min;
        bool is_3d_flag;
    public:
        Accuracy_Atm(int im, int jm, int km, double dthe, double dphi);
        Accuracy_Atm(int im, int jm, int km, double dr, double dthe, double dphi);
        ~Accuracy_Atm();
        std::tuple<double, int, int> residuumQuery_2D( 
            Array_1D &rad, Array_1D &the, Array &v, Array &w, Vector3D<> & residuum_2d);
        std::tuple<double, int, int, int> 
            residuumQuery_3D(Array_1D &rad, Array_1D &the, Array &u, Array &v,
            Array &w, Vector3D<> & residuum_3d);
        void steadyQuery_2D(Array &v, Array &vn, Array &w, Array &wn,
            Array &p_dyn, Array &p_dynn);
        void steadyQuery_3D(Array &u, Array &un, Array &v, Array &vn,
            Array &w, Array &wn, Array &t, Array &tn, Array &c, Array &cn,
            Array &cloud, Array &cloudn, Array &ice, Array &icen, Array &co2,
            Array &co2n, Array &p_dyn, Array &p_dynn, double L_atm);
        void print(const string& name, double value, int j, int k) const;        
        void print(const string& name, double value, int i, int j, int k) const;
        double out_min() const;
        int out_i_res() const;
        int out_j_res() const;
        int out_k_res() const;
};
#endif

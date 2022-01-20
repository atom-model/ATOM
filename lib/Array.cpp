#include <algorithm>
#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>

#include "Array.h"
#include "Utils.h"
#include "cAtmosphereModel.h"
#include "cHydrosphereModel.h"

using namespace std;
using namespace AtomUtils;

#define MAXI 41
#define MAXJ 181
#define MAXK 361

Array::Array(int idim, int jdim, int kdim, double val):
    x(NULL){
        initArray(idim, jdim, kdim, val);
}
/*
*
*/
Array::~Array(){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            delete[] x[ i ][ j ];
        }
    }
    for(int i = 0; i < im; i++){
        delete[] x[ i ];
    }
    delete[] x;
}
/*
*
*/
void Array::initArray(int im, int jm, int km, double aa){
    if(x){
        assert(im == this->im);
        assert(jm == this->jm);
        assert(km == this->km);
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                for(int k = 0; k < km; k++){
                    x[ i ][ j ][ k ] = aa;
                }
            }
        }
    }else{
        assert(im <= MAXI);
        assert(jm <= MAXJ);
        assert(km <= MAXK);
        this->im = im;
        this->jm = jm;
        this->km = km;
        x = new double**[im];
        for(int i = 0; i < im; i++){
            x[ i ] = new double*[jm];
            for(int j = 0; j < jm; j++){
                x[ i ][ j ] = new double[km];
                for(int k = 0; k < km; k++){
                    x[ i ][ j ][ k ] = aa;
                }
            }
        }    
    }
}
/*
*
*/
void Array::printArray(string dir, int im, int jm, int km){
    assert(im == this->im);
    assert(jm == this->jm);
    assert(km == this->km);

    cout.precision(6);
    cout.setf(ios::fixed);

    int i_beg = 0;
    int i_end = 0;
    string at = "AGCM";
    string hy = "OGCM";
    if(dir.compare(at) == 0){
        i_beg = 0;
        i_end = 0;
        cout << endl << endl << " ...................... \
 printout from atmosphere code at level = " << i_beg << "  ...................." 
            << endl << endl;
    }
    if(dir.compare(hy) == 0){
        i_beg = im-1;
        i_end = im-1;
        cout << endl << endl << " ...................... \
 printout from hydrosphere code at level = " << i_end << "  ...................." 
            << endl << endl;
    }
    for(int i = i_beg; i <= i_end; i++){
        cout << endl;
        cout << "  phi = k-direction in  " << km-1 
            << "  longitudinal steps " << endl
             << "  theta = j-direction in  " << jm-1 
             << "  latitudinal steps " << endl
             << "  r = i-direction at level  " << i << endl << endl;
        for(int j = 0; j < jm; j+=4){
            for(int k = 0; k < km; k+=20){
                cout.width(4);
                cout.fill(' ');
                cout << x[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}
/*
*
*/
void Array::inspect(const std::string& prefix) const{
    std::vector<double> mins(im, 0), maxes(im, 0), means(im, 0), s_means(im, 0);
    for(int i=0; i<im; i++){
        double min_tmp=x[i][0][0], max_tmp=x[i][0][0], mean_tmp=0, s_means_tmp=0, weight_tmp=0;
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                min_tmp=std::min(min_tmp, x[i][j][k]);
                max_tmp=std::max(max_tmp, x[i][j][k]);
                mean_tmp+=x[i][j][k];
                double w=cos(abs(90-j)*M_PI/180.);
                s_means_tmp+=x[i][j][k]*w;
                weight_tmp+=w;
            }
        }
        mins[i]=min_tmp;
        maxes[i]=max_tmp;
        means[i]=mean_tmp/(jm*km);
        s_means[i]=s_means_tmp/weight_tmp;
    }
    logger()<<prefix<<"==================================="<<std::endl;
    logger()<<prefix<<"max:: " << *std::max_element(maxes.begin(), maxes.end()) << std::endl;
    logger()<<prefix<< " ";
    for(std::vector<double>::iterator it=maxes.begin(); it!=maxes.end(); it++){
        logger()<< fixed << setprecision(4) << *it << "  ";
    }
    logger()<<std::endl;
    logger()<<prefix<<"min:: " << *std::min_element(mins.begin(), mins.end()) << std::endl;
    logger()<<prefix<< " ";
    for(std::vector<double>::iterator it=mins.begin(); it!=mins.end(); it++){
        logger()<< *it << "  ";
    }
    logger()<<std::endl;
    logger()<<prefix<<"mean:: " << std::accumulate(means.begin(), means.end(), 0.0)/means.size() << std::endl;
    logger()<<prefix<< " ";
    for(std::vector<double>::iterator it=means.begin(); it!=means.end(); it++){
        logger()<< *it << "  ";
    }
    logger()<<std::endl;
    logger()<<prefix<<"spherical mean of each layer:: " << std::endl;
    logger()<<prefix<< " ";
    for(std::vector<double>::iterator it=s_means.begin(); it!=s_means.end(); it++){
        logger()<< *it << "  ";
    }
    logger()<<std::endl;
    logger()<<prefix<<"==================================="<<std::endl;
}


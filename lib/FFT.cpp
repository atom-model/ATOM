#include <complex>
#include <iostream>
#include <iomanip>
#include <valarray>
#include <fstream>
#include <math.h>       

#include <Array.h>
#include <Utils.h>

//#define FFT_DEBUG

//Some of the code below was copied from
//https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B
 
const double PI = 3.141592653589793238460;
 
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
 
// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];
 
    // conquer
    fft(even);
    fft(odd);
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}
 
// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
// !!! Warning : in some cases this code make result different from not optimased version above (need to fix bug)
// The bug is now fixed @2017/05/30 
void fft_ex(CArray &x)
{
    // DFT
    unsigned int N = x.size(), k = N, n;
    double thetaT = 3.14159265358979323846264338328L / N;
    Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
    while (k > 1)
    {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0L;
        for (unsigned int l = 0; l < k; l++)
        {
            for (unsigned int a = l; a < N; a += n)
            {
                unsigned int b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    // Decimate
    unsigned int m = (unsigned int)log2(N);
    for (unsigned int a = 0; a < N; a++)
    {
        unsigned int b = a;
        // Reverse bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a)
        {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }
    //// Normalize (This section make it not working correctly)
    //Complex f = 1.0 / sqrt(N);
    //for (unsigned int i = 0; i < N; i++)
    //  x[i] *= f;
}
 
// inverse fft (in-place)
void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);
 
    // forward fft
    fft_ex( x );
 
    // conjugate the complex numbers again
    x = x.apply(std::conj);
 
    // scale the numbers
    x /= x.size();
}

int next_power_of_two(unsigned number)
{
    unsigned mask = 0x80000000;
    for (int i = 31; i >= 0; i--)
    {
        if (number & mask) 
            return 1 << (i+1);
        else
            mask >>= 1;
    }
    return 0;
}

void init_gaussian_kernel(double* kernel, double sigma, size_t len)
{
    double t = sqrt(2*M_PI)*sigma, ss=2*sigma*sigma;
    //std::cout << t << " " << ss <<std::endl;
    double sum=0;
    for(size_t i=0; i<len; i++)
    {
        //std::cout << i << std::endl;
        double x_2 = (i-len/2)*(i-len/2);
        kernel[i] = exp(-x_2/ss)/t;
        sum += kernel[i]; 
        //std::cout << kernel[i] << std::endl;
    }
    for(size_t i=0; i<len; i++)
    { 
        kernel[i] /= sum; 
    }
}

void AtomUtils::fft_gaussian_filter(Array& in_data, int sigma){
    using namespace std;
    int jm=in_data.jm, km=in_data.km, im=in_data.im;

    int nn = next_power_of_two(jm);
    double* kernel = new double[nn];
    init_gaussian_kernel(kernel, sigma, nn);
    CArray k_data(nn);
    for(int i =0; i<nn; i++){
        //std::cout << kernel[i] << ", ";
        k_data[i] = Complex(kernel[i]);
    }
    //std::cout << std::endl;
    delete[] kernel;

    //fft gaussian kernel data
    fft_ex(k_data);

    for(int i=0; i<im; i++){
        for(int k=0; k<km; k++)
        {
            CArray data(nn);
            for(int j=0; j<nn; j++)
            {
                if(j < jm)
                    data[j] = Complex(in_data.x[i][j][k]);
                else{//mirror the padding
                    if(j<jm+(nn-jm)/2)
                        data[j] = Complex(in_data.x[i][jm-(j-jm)][k]);
                    else
                        data[j] = Complex(in_data.x[i][nn-j][k]);
                }
            }
 
            // forward fft
            fft_ex(data);

            //apply the gaussian filter in frequency domain
            for(int j =0; j<nn; j++)
            {
                Complex tmp(data[j].real()*k_data[j].real(), data[j].imag()*k_data[j].real());
                data[j] = tmp;
            }

            // inverse fft
            ifft(data);

            for(int j =0; j<nn; j++)
            {
                //the result data need shift
                if(j>=nn/2){
                    in_data.x[i][j-nn/2][k] = data[j].real();
                }else{
                    if(j+nn/2<jm)
                        in_data.x[i][j+nn/2][k] = data[j].real();
                }
            } 
        }  
    }
    return;
}

#ifdef FFT_DEBUG
 
int main()
{
    using namespace std;
    int m=361, n=181;

    int nn = next_power_of_two(n);
    double* kernel = new double[nn];
    init_gaussian_kernel(kernel, 5, nn);
    CArray k_data(nn);
    for(int i =0; i<nn; i++){
        std::cout << kernel[i] << ", ";
        k_data[i] = Complex(kernel[i]);
    }
    std::cout << std::endl;
    delete[] kernel;

    fft_ex(k_data);

    /*std::cout << "kernel data" << std::endl;
    for (int i = 0; i < n; ++i)
    {
        std::cout << k_data[i] << " ";
    }
    std::cout << std::endl << std::endl;
    */
    double* src = new double[m*n](); 
    double* dst = new double[m*nn]();
    std::fill_n(src+15*m, m*15, 1);
    std::fill_n(src+45*m, m*15, 0.2);
    std::fill_n(src+75*m, m*15, 0.4);
    std::fill_n(src+105*m, m*15, 0.6);
    std::fill_n(src+135*m, m*15, 0.8);
    std::fill_n(src+160*m, m*10, 1);

    std::ofstream ios("input.bin", std::ios::binary | std::ios::out);
    for(int j=0; j<m*n; j++){
        ios.write(reinterpret_cast<char*>(&src[j]), sizeof(double)) ;
    }
    ios.close();

    for(int j=0; j<m; j++)
    { 
        CArray data(nn);
        for(int i =0; i<nn; i++)
        {
            if(i < n)
                data[i] = Complex(src[m*i+j]);
            else
                data[i] = Complex(0);
        }

        std::cout << std::fixed;
        std::cout << std::setprecision(2);
    
        /*std::cout << "input data" << std::endl;
        for (int i = 0; i < 8; ++i)
        {
            std::cout << data[i] << std::endl;
        }*/

        // forward fft
        fft_ex(data);

        /*std::cout << std::endl << "fft" << std::endl;
        for (int i = 0; i < 8; ++i)
        {
            std::cout << data[i] << std::endl;
        }*/

        for(int i =0; i<nn; i++)
        {
            //double x=data[i].real(), y=data[i].imag(), u=k_data[i].real(), v=k_data[i].imag();
            //Complex tmp=(x*u-y*v, x*v+y*u);
            Complex tmp(data[i].real()*k_data[i].real(), data[i].imag()*k_data[i].real());
            data[i] = tmp;
        }
 
        // inverse fft
        ifft(data);
 
        /*std::cout << std::endl << "ifft" << std::endl;
        for (int i = 0; i < 8; ++i)
        {
            std::cout << data[i] << std::endl;
        }*/

        for(int i =0; i<nn; i++)
        {
            //the result data need shift
            if(i>=nn/2){
                dst[m*(i-nn/2)+j] = data[i].real();
            }else{
                dst[m*(i+nn/2)+j] = data[i].real();
            }
            //dst[m*i+j] = data[i].real();
            if(j==0)
                std::cout << data[i].real() << ",";
        }
    }

    std::ofstream os("output.bin", std::ios::binary | std::ios::out);
    for(int j=0; j<m*n; j++){
        os.write(reinterpret_cast<char*>(&dst[j]), sizeof(double)) ;
    }
    os.close();

    delete[] src;
    delete[] dst;

    return 0;
}
#endif

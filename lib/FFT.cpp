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

    assert(nn<2*jm);

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
                        data[j] = Complex(in_data.x[i][(jm-1)-(j-jm)][k]);
                    else
                        data[j] = Complex(in_data.x[i][(nn-1)-j][k]);
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

//mirror padding the data
//i_len: the length of the input data
//p_len: the length after the padding
//data:  array of size p_len, the first i_len are input data, this function will pad it to p_len
void AtomUtils::mirror_padding(double* data, size_t i_len, size_t p_len)
{
    if (p_len <= i_len) return;

    for(size_t j=i_len; j<p_len; j++)
    {
        if(j<i_len+(p_len-i_len)/2)//the first half
            data[j] = data[(i_len-1)-(j-i_len)]; //use the data from the end
        else//second half
            data[j] = data[(p_len-1)-j];//use the data from the beginning
    }
}

//data and kernel must be the same size len
//result will replace the input data
void AtomUtils::fft_gaussian_filter(double* _data, double* kernel, size_t len){
    CArray k_data(len), data(len);
    for(size_t i =0; i<len; i++){
        k_data[i] = Complex(kernel[i]);
        data[i] = Complex(_data[i]);
    }

    //fft gaussian kernel data
    fft_ex(k_data);
    //fft input data
    fft_ex(data);

    //apply the gaussian filter in frequency domain
    for(size_t j =0; j<len; j++)
    {
        data[j] = Complex(data[j].real()*k_data[j].real(), data[j].imag()*k_data[j].real());
    }

    // inverse fft
    ifft(data);

    for(size_t j =0; j<len; j++)
    {
        //the result data need shift
        if(j>=len/2){
            _data[j-len/2] = data[j].real();
        }else{
            _data[j+len/2] = data[j].real();
        }
    }
    return; 
}

//do the fft in 3 directions
void AtomUtils::fft_gaussian_filter_3d(Array& data, int sigma)
{
    using namespace std;
    int jm=data.jm, km=data.km, im=data.im;

    //filter k direction
    int len_k = next_power_of_two(km); //data size must be power of 2
    //init kernel
    double* kernel_k = new double[len_k];
    init_gaussian_kernel(kernel_k, sigma, len_k);

    double* src_k = new double[len_k]();
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            //take the data out    
            for(int k=0; k<km; k++){
                src_k[k] = data.x[i][j][k]; 
            }
            //padding and filter
            mirror_padding(src_k, km, len_k);
            fft_gaussian_filter(src_k, kernel_k, len_k);
            //put the filtered data back
            for(int k=0; k<km; k++){
                data.x[i][j][k] = src_k[k]; 
            }
        }
    }
    delete[] kernel_k;
    delete[] src_k;


    //filter i direction
    int len_i = next_power_of_two(im); //data size must be power of 2
    //init kernel
    double* kernel_i = new double[len_i];
    init_gaussian_kernel(kernel_i, sigma, len_i);

    double* src_i = new double[len_i]();
    for(int k=0; k<km; k++){
        for(int j=0; j<jm; j++){
            //take the data out    
            for(int i=0; i<im; i++){
                src_i[i] = data.x[i][j][k];
            }
            //padding and filter
            mirror_padding(src_i, im, len_i);
            fft_gaussian_filter(src_i, kernel_i, len_i);
            //put the filtered data back
            for(int i=0; i<im; i++){
                data.x[i][j][k] = src_i[i];
            }
        }
    }
    delete[] kernel_i;
    delete[] src_i;


    //filter j direction
    int len_j = next_power_of_two(jm); //data size must be power of 2
    //init kernel
    double* kernel_j = new double[len_j];
    init_gaussian_kernel(kernel_j, sigma, len_j);

    double* src_j = new double[len_j]();
    for(int i=0; i<im; i++){
        for(int k=0; k<km; k++){
            //take the data out    
            for(int j=0; j<jm; j++){
                src_j[j] = data.x[i][j][k];
            }
            //padding and filter
            mirror_padding(src_j, jm, len_j);
            fft_gaussian_filter(src_j, kernel_j, len_j);
            //put the filtered data back
            for(int j=0; j<jm; j++){
                data.x[i][j][k] = src_j[j];
            }
        }
    }
    delete[] kernel_j;
    delete[] src_j;

}

#ifdef FFT_DEBUG
 
int main()
{
    using namespace std;
    using namespace AtomUtils;
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
    std::cout << std::endl;
    std::ofstream os("output.bin", std::ios::binary | std::ios::out);
    for(int j=0; j<m*n; j++){
        os.write(reinterpret_cast<char*>(&dst[j]), sizeof(double)) ;
    }
    os.close();

    delete[] src;
    delete[] dst;

    //test fft 1-d array
    int len_1 = next_power_of_two(100);
    double* kernel_1 = new double[len_1];
    
    init_gaussian_kernel(kernel_1, 5, len_1);
    
    double* src_1 = new double[len_1](); 
    for(int i=0; i<10; i++)
        std::fill_n(src_1+10*i, 10, i);

    std::cout << "original data: " << std::endl;
    for(int i=0; i<len_1; i++)
        std::cout << src_1[i] << ", ";
    std::cout << std::endl;

    mirror_padding(src_1, 100, len_1);
    fft_gaussian_filter(src_1, kernel_1, len_1);

    std::cout << "filtered data: " << std::endl;
    for(int i=0; i<len_1; i++)
        std::cout << src_1[i] << ", ";
    std::cout << std::endl;

    delete[] kernel_1;
    delete[] src_1;

    return 0;
}
#endif

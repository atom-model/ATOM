#include <complex>
#include <iostream>
#include <iomanip>
#include <valarray>
#include <fstream>
#include <math.h>       
 
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
    fft( x );
 
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
 
int main()
{
    using namespace std;

    std::ifstream input_file ("./atm_t_time_0_iter_n_layer_0.bin", ios::in | ios::binary);
    input_file.seekg(0, std::ios::end);
    size_t size = input_file.tellg();  
    input_file.seekg(0, std::ios::beg); 
    char* buffer = new char[size];
    input_file.read(buffer, size);
    input_file.close();

    size_t size_n = next_power_of_two(size/8);
    std::cout << size/8 << " " << log2(size_n) << " " << size_n << std::endl;

    double* padded_double_data = new double[size_n](); 
    std::copy((double*)buffer, (double*)buffer+size/8, padded_double_data);
    delete[] buffer;
    //const Complex test[] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
    CArray data(size_n);
    for(int i =0; i<size_n; i++){
        data[i] = Complex(padded_double_data[i]);
    }
    delete[] padded_double_data;
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    
    std::cout << "input data" << std::endl;
    for (int i = 0; i < 8; ++i)
    {
        std::cout << data[i] << std::endl;
    }

    // forward fft
    fft_ex(data);

    std::cout << std::endl << "fft" << std::endl;
    for (int i = 0; i < 8; ++i)
    {
        std::cout << data[i] << std::endl;
    }
 
    // inverse fft
    ifft(data);
 
    std::cout << std::endl << "ifft" << std::endl;
    for (int i = 0; i < 8; ++i)
    {
        std::cout << data[i] << std::endl;
    }

    std::ofstream os("test.bin", std::ios::binary | std::ios::out);
    for(int j=0; j<size/8; j++){
        os.write(reinterpret_cast<char*>(&data[j].real()), sizeof(double)) ;
    }
    os.close();
    return 0;
}

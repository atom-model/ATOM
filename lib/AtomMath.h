/*#include <iostream>

double parabola_interp(double lower_bound, double upper_bound, double x);

int main(){
    std::cout << "hello world" << std::endl;
    std::cout << parabola_interp(10, 20, 0) << std::endl;
    std::cout << parabola_interp(10, 20, 0.25) << std::endl;
    std::cout << parabola_interp(10, 20, 0.5) << std::endl;
    std::cout << parabola_interp(10, 20, 0.75) << std::endl;
    std::cout << parabola_interp(10, 20, 1) << std::endl;
    std::cout << parabola_interp(10, 20, 1.25) << std::endl;
    std::cout << parabola_interp(10, 20, 1.5) << std::endl;
    std::cout << parabola_interp(10, 20, 1.75) << std::endl;
    std::cout << parabola_interp(10, 20, 2) << std::endl;
}
*/
//x^2-2x: x[0,2], y[-1,0]
//lower_bound: the value when x=1(y=-1)
//upper_bound: the value when x=0 or 2(y=0)
//for example parabola_interp(10, 20, 1) return 10
//parabola_interp(10, 20, 0) return 20
double parabola_interp(double lower_bound, double upper_bound, double x){
    return upper_bound + (upper_bound-lower_bound)*(x*x - 2*x);
}


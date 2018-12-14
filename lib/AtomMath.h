#ifndef _ATOMMATH_
#define _ATOMMATH_

//x^2-2x: x[0,2], y[-1,0]
//lower_bound: the value when x=1(y=-1)
//upper_bound: the value when x=0 or 2(y=0)
//for example parabola_interp(10, 20, 1) return 10
//parabola_interp(10, 20, 0) return 20
inline double parabola_interp(double lower_bound, double upper_bound, double x){
    return upper_bound + (upper_bound-lower_bound)*(x*x - 2*x);
}

#endif

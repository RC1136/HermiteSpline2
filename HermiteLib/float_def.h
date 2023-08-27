
#ifndef FLOAT_DEF_H
#define FLOAT_DEF_H

#include <math.h>

// #define __LONG_DOUBLE
// #define __FLOAT_128

#if defined __LONG_DOUBLE


typedef long double myfloat_t;

#define acos(x)    acosl(x)
#define asin(x)    asinl(x)
#define atan(x)    atanl(x)
#define atan2(y,x) atan2l(x)
#define cos(x)     cosl(x)
#define cosh(x)    coshl(x)
#define sin(x)     sinl(x)
#define sinh(x)    sinhl(x)
#define tanh(x)    tanhl(x)
#define frexp(x,e) frexpl(x,e)
#define ldexp(x,e) ldexpl(x,e)
#define log(x)     logl(x)
#define log10(x)   log10l(x)
#define modf(x,i)  modfl(x,i)
#define pow(x,y)   powl(x,y)
#define sqrt(x)    sqrtl(x)
#define ceil(x)    ceill(x)
#define fabs(x)    fabsl(x)
#define floor(x)   floorl(x)
#define fmod(x,y)  fmodl(x,y)
#define fmax(x,y)  fmaxl(x,y)
#define fmin(x,y)  fminl(x,y)

#elif defined __FLOAT_128

#include <quadmath.h>

typedef __float128 myfloat_t;

#define acos(x)    acosq(x)
#define asin(x)    asinq(x)
#define atan(x)    atanq(x)
#define atan2(y,x) atan2q(x)
#define cos(x)     cosq(x)
#define cosh(x)    coshq(x)
#define sin(x)     sinq(x)
#define sinh(x)    sinhq(x)
#define tanh(x)    tanhq(x)
#define frexp(x,e) frexpq(x,e)
#define ldexp(x,e) ldexpq(x,e)
#define log(x)     logq(x)
#define log10(x)   log10q(x)
#define modf(x,i)  modfq(x,i)
#define pow(x,y)   powq(x,y)
#define sqrt(x)    sqrtq(x)
#define ceil(x)    ceilq(x)
#define fabs(x)    fabsq(x)
#define floor(x)   floorq(x)
#define fmod(x,y)  fmodq(x,y)
#define fmax(x,y)  fmaxq(x,y)
#define fmin(x,y)  fminq(x,y)

#undef isnan
#define isnan(x)   isnanq((myfloat_t)(x))
#undef isinf
#define isinf(x)   isinfq((myfloat_t)(x))

#else

typedef double myfloat_t;

#endif


#endif

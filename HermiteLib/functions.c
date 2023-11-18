#include <math.h>
#include "functions.h"
#include "hermite.h"

myfloat_t testfunc0(const myfloat_t _x)
{
	//a = 1, b = -3, c = 4, d = 7
	//return 0.0;
	return W21(_x, (myfloat_t[]) { 1, -3, 4, 7 });
}

myfloat_t testdf0(const myfloat_t _x)
{
	//return 0.0;
	return dW21(_x, (myfloat_t[]) { 1, -3, 4, 7 });
}

myfloat_t testfunc1(const myfloat_t _x)
{
	//a = 1, b = -2, c = 5, d = 3, v = 2
	//return 0.0;
	return W31(_x, (myfloat_t[]) { 1, -2, 5, 3, 2 });
}

myfloat_t testdf1(const myfloat_t _x)
{
	//return 0.0;
	return dW31(_x, (myfloat_t[]) { 1, -2, 5, 3, 2 });
}

myfloat_t testfunc2(const myfloat_t _x)
{
	//a = 3, b = 2, c = -4, d = 1
	//return 0.0;
	return 3 + 2 * _x - 4 * _x * _x + 1 * _x * _x * _x;
}

myfloat_t testdf2(const myfloat_t _x)
{
	return 2 - 8 * _x + 3 * _x * _x;
}

myfloat_t testfunc3(const myfloat_t _x)
{
	//a = 3, b = 2, c = -4, d = 1, h = -1
	return 3 + 2 * _x - 4 * _x * _x + 1 * _x * _x * _x - 1 * _x * _x * _x * _x;
}

myfloat_t testdf3(const myfloat_t _x)
{
	return 2 - 8 * _x + 3 * _x * _x - 4 * _x * _x * _x;
}

myfloat_t testfunc4(const myfloat_t _x)
{
	return exp(pow(_x - 4, 3) / 42);
}

myfloat_t testdf4(const myfloat_t _x)
{
	return exp(pow(_x - 4., 3.) / 42.) * pow(_x - 4., 2.) / 14.;
}

myfloat_t testfunc5(const myfloat_t _x)
{
	return sin(_x);
}

myfloat_t testdf5(const myfloat_t _x)
{
	return cos(_x);
}

myfloat_t testfunc6(const myfloat_t _x)
{
	return sqrt(pow(_x, 3.) + 1.);
}

myfloat_t testdf6(const myfloat_t _x)
{
	return (3. / 2.) * ((_x * _x) / sqrt(_x * _x * _x + 1));
}

myfloat_t testfunc7(const myfloat_t _x)
{
	return 1. / (_x * _x + 0.4);
}

myfloat_t testdf7(const myfloat_t _x)
{
	return (-2 * _x) / pow(_x * _x + 0.4, 2);
}

myfloat_t testfunc8(const myfloat_t _x)
{
	return tan(_x) / (_x + 2.) + 2;
}

myfloat_t testdf8(const myfloat_t _x)
{
	return (1 + pow(tan(_x), 2)) / (_x + 2) - (tan(_x)) / pow(_x + 2, 2);
}

myfloat_t testfunc9(const myfloat_t _x)
{
	return log(_x);
}

myfloat_t testdf9(const myfloat_t _x)
{
	return 1 / _x;
}

myfloat_t testfunc10(const myfloat_t _x)
{
	return exp(sin(_x) - cos(_x));
}

myfloat_t testdf10(const myfloat_t _x)
{
	return (cos(_x) + sin(_x)) * testfunc10(_x);
}

myfloat_t testfunc11(const myfloat_t _x)
{
	return exp(sin(_x)) + _x;
}

myfloat_t testdf11(const myfloat_t _x)
{
	return cos(_x) * exp(sin(_x)) + 1;
}

myfloat_t testfunc12(const myfloat_t _x)
{
	return tan(_x);
}

myfloat_t testdf12(const myfloat_t _x)
{
	return 1. / pow(cos(_x), 2);
}

myfloat_t testfunc13(const myfloat_t _x)
{
	return 1. / (_x * _x * _x + 1);
}

myfloat_t testdf13(const myfloat_t _x)
{
	return (-3 * _x * _x) / pow(_x * _x * _x + 1, 2);
}

myfloat_t testfunc14(const myfloat_t _x)
{
	return 1. / (1.0 + _x * _x);
}

myfloat_t testdf14(const myfloat_t _x)
{
	return (-2 * _x) / pow(_x * _x + 1., 2);
}

myfloat_t testfunc15(const myfloat_t _x)
{
	return cos(_x);
}

myfloat_t testdf15(const myfloat_t _x)
{
	return -sin(_x);
}

myfloat_t testfunc16(const myfloat_t _x)
{
	return exp(sin(_x));
}

myfloat_t testdf16(const myfloat_t _x)
{
	return exp(sin(_x)) * cos(_x);
}

myfloat_t testfunc17(const myfloat_t _x)
{
	return asin(_x);
}

myfloat_t testdf17(const myfloat_t _x)
{
	return 1. / sqrt(1 - _x*_x);
}

myfloat_t testfunc18(const myfloat_t _x)
{
	return sinh(_x);
}

myfloat_t testdf18(const myfloat_t _x)
{
	return cosh(_x);
}

myfloat_t testfunc19(const myfloat_t _x)
{
	return cosh(_x);
}

myfloat_t testdf19(const myfloat_t _x)
{
	return sinh(_x);
}


const function funcs[] =  { testfunc0, testfunc1, testfunc2, testfunc3, testfunc4, testfunc5, testfunc6, testfunc7, testfunc8, testfunc9, testfunc10, testfunc11, testfunc12, testfunc13, testfunc14, testfunc15, testfunc16, testfunc17, testfunc18, testfunc19, };
const function dfuncs[] = { testdf0,   testdf1,   testdf2,   testdf3,   testdf4,   testdf5,   testdf6,   testdf7,   testdf8,   testdf9,   testdf10,   testdf11,   testdf12,   testdf13,   testdf14,   testdf15,   testdf16,   testdf17,   testdf18,   testdf19, };

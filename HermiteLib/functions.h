
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "float_def.h"

typedef myfloat_t (*function)(const myfloat_t);


extern myfloat_t testfunc0(const myfloat_t);
extern myfloat_t testfunc1(const myfloat_t);
extern myfloat_t testfunc2(const myfloat_t);
extern myfloat_t testfunc3(const myfloat_t);
extern myfloat_t testfunc4(const myfloat_t);
extern myfloat_t testfunc5(const myfloat_t);
extern myfloat_t testfunc6(const myfloat_t);
extern myfloat_t testfunc7(const myfloat_t);
extern myfloat_t testfunc8(const myfloat_t);
extern myfloat_t testfunc9(const myfloat_t);
extern myfloat_t testfunc10(const myfloat_t);
extern myfloat_t testfunc11(const myfloat_t);
extern myfloat_t testfunc12(const myfloat_t);
extern myfloat_t testfunc13(const myfloat_t);
extern myfloat_t testdf0(const myfloat_t);
extern myfloat_t testdf1(const myfloat_t);
extern myfloat_t testdf2(const myfloat_t);
extern myfloat_t testdf3(const myfloat_t);
extern myfloat_t testdf4(const myfloat_t);
extern myfloat_t testdf5(const myfloat_t);
extern myfloat_t testdf6(const myfloat_t);
extern myfloat_t testdf7(const myfloat_t);
extern myfloat_t testdf8(const myfloat_t);
extern myfloat_t testdf9(const myfloat_t);
extern myfloat_t testdf10(const myfloat_t);
extern myfloat_t testdf11(const myfloat_t);
extern myfloat_t testdf12(const myfloat_t);
extern myfloat_t testdf13(const myfloat_t);




extern const function funcs[];
extern const function dfuncs[];

#endif

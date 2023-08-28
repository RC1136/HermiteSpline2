
#ifndef HERMITE_H
#define HERMITE_H

#include "float_def.h"

//Вигляд ланки
enum linktype {
	powexp4  = 1,
	powexp5  = 2,
	poly4    = 3,
	poly5    = 4,
	exppow5  = 5,
	pow1exp2 = 6,
	pow2exp2 = 7,
	linktype_count
};

typedef struct {
	enum linktype type;
	int param_count; //к-сть параметрів у одній ланці
	int link_count;	 //кількість ланок
	double* A;		 //параметри сплайна (link_conut*param_count елементів)
	double* X;		 //точки наближення (link_count+1 елементів)
	myfloat_t* A128;
	myfloat_t* X128;
} herm_params;

typedef myfloat_t(*function)(const myfloat_t);

myfloat_t PE4(const myfloat_t x, const myfloat_t a[4]);

myfloat_t dPE4(const myfloat_t x, const myfloat_t a[4]);

myfloat_t PE5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dPE5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t EP5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dEP5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t PN4(const myfloat_t x, const myfloat_t a[4]);

myfloat_t dPN4(const myfloat_t x, const myfloat_t a[4]);

myfloat_t PN5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dPN5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t W22(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dW22(const myfloat_t x, const myfloat_t a[5]);

myfloat_t W12(const myfloat_t x, const myfloat_t a[4]);

myfloat_t dW12(const myfloat_t x, const myfloat_t a[4]);

myfloat_t HermiteSpline(const herm_params hp, const myfloat_t x, const char derivative);

int HermGen(function _f[], herm_params* hp, const myfloat_t a, const myfloat_t b, const myfloat_t nu);

int HermGen2(function _f[], herm_params* hp, const myfloat_t a, const myfloat_t b, const int r);

myfloat_t finderr(myfloat_t (*link)(const myfloat_t, const myfloat_t[]), const myfloat_t params[], function f, const myfloat_t from, const myfloat_t to);

#endif

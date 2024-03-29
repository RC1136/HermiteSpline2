
#ifndef HERMITE_H
#define HERMITE_H

#include "float_def.h"

//������ �����
enum linktype {
	pow2exp1  = 1,
	pow3exp1  = 2,
	poly4    = 3,
	poly5    = 4,
	pow1exp3  = 5,
	pow1exp2 = 6,
	pow2exp2 = 7,
	linktype_count
};

typedef struct {
	enum linktype type;
	int param_count; //�-��� ��������� � ����� �����
	int link_count;	 //������� �����
	double* A;		 //��������� ������� (link_conut*param_count ��������)
	double* X;		 //����� ���������� (link_count+1 ��������)
	myfloat_t* A128;
	myfloat_t* X128;
} herm_params;

typedef myfloat_t(*function)(const myfloat_t);

myfloat_t W21(const myfloat_t x, const myfloat_t a[4]);

myfloat_t dW21(const myfloat_t x, const myfloat_t a[4]);

myfloat_t W31(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dW31(const myfloat_t x, const myfloat_t a[5]);

myfloat_t W13(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dW13(const myfloat_t x, const myfloat_t a[5]);

myfloat_t PN4(const myfloat_t x, const myfloat_t a[4]);

myfloat_t dPN4(const myfloat_t x, const myfloat_t a[4]);

myfloat_t PN5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dPN5(const myfloat_t x, const myfloat_t a[5]);

myfloat_t W22(const myfloat_t x, const myfloat_t a[5]);

myfloat_t dW22(const myfloat_t x, const myfloat_t a[5]);

myfloat_t W12(const myfloat_t x, const myfloat_t a[4]);

myfloat_t dW12(const myfloat_t x, const myfloat_t a[4]);

myfloat_t HermiteSpline(const herm_params hp, const myfloat_t x, const unsigned char derivative);

int HermGenNu(function _f[], herm_params* hp, const myfloat_t a, const myfloat_t b, const myfloat_t nu);

int HermGenR(function _f[], herm_params* hp, const myfloat_t a, const myfloat_t b, const int r);

myfloat_t finderr(myfloat_t (*link)(const myfloat_t, const myfloat_t[]), const myfloat_t params[], function f, const myfloat_t from, const myfloat_t to);

void HermClear(herm_params* hp);

#endif

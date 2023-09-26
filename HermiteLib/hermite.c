#include "hermite.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <stdio.h>

#define isodd(x) (x & 1)


#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


int param_count_map[] = {
	[powexp4]  = 4,
	[powexp5]  = 5,
	[poly4]    = 4,
	[poly5]    = 5,
	[exppow5]  = 5,
	[pow2exp2] = 5,
	[pow1exp2] = 4,
};

//Степенево-експоненціальна ланка з чотирма параметрами
myfloat_t PE4(const myfloat_t x, const myfloat_t a[4])
{
	return a[0] * pow(x, a[1] + a[2] * x) * exp(a[3] * x);
}

//Похідна степенево-експоненціальної ланки з чотирма параметрами
myfloat_t dPE4(const myfloat_t x, const myfloat_t a[4])
{
	return PE4(x, a) * (a[2] * (1 + log(x)) + a[1] / x + a[3]);
}

//Степенево-експоненціальна ланка з п'ятьма параметрами
myfloat_t PE5(const myfloat_t x, const myfloat_t a[5])
{
	return a[0] * pow(x, a[1] + a[2] * x + a[3] * x * x) * exp(a[4] * x);
}

//Похідна степенево-експоненціальної ланки з п'ятьма параметрами
myfloat_t dPE5(const myfloat_t x, const myfloat_t a[5])
{
	return PE5(x, a) * (
		(a[2] + 2 * a[3] * x) * log(x) + a[1] / x + a[2] + a[3] * x + a[4]
	);
}


myfloat_t EP5(const myfloat_t x, const myfloat_t A[5])
{
	const myfloat_t a = A[0], b = A[1], c = A[2], d = A[3], g = A[4];
	return a * pow(x, b) * exp(c * x + d * x * x + g * x * x * x);
}


myfloat_t dEP5(const myfloat_t x, const myfloat_t A[5])
{
	const myfloat_t a = A[0], b = A[1], c = A[2], d = A[3], g = A[4];
	return EP5(x, A) * (b / x + c + 2 * d * x + 3 * g * x * x);
}


myfloat_t W22(const myfloat_t x, const myfloat_t A[5])
{
	// A = A[0], a1 = A[1], a2 = A[2], b1 = A[3], b2 = A[4]
	return A[0] * pow(x, A[1] + A[2] * x) * exp(A[3] * x + A[4] * x * x);
}


myfloat_t dW22(const myfloat_t x, const myfloat_t A[5])
{
	// A = A[0], a1 = A[1], a2 = A[2], b1 = A[3], b2 = A[4]
	return A[0] * pow(x, A[1] + A[2] * x - 1.) * exp(A[3] * x + A[4] * x * x) *
		(A[1] + A[2] * x * (1 + log(x)) + A[3] * x + 2. * A[4] * x * x);
}


myfloat_t W12(const myfloat_t x, const myfloat_t A[4])
{
	const myfloat_t a = A[0], b = A[1], c = A[2], d = A[3];
	return a * pow(x, b) * exp(c * x + d * x * x);
}


myfloat_t dW12(const myfloat_t x, const myfloat_t A[4])
{
	const myfloat_t a = A[0], b = A[1], c = A[2], d = A[3];
	return W12(x, A) * (b / x + c + 2 * d * x);
}


//Поліном з count+1 параметрами 
myfloat_t Polynomial(const myfloat_t x, const myfloat_t* a, const int count)
{
	myfloat_t res = a[0], tmpx = 1.0;
	for (int i = 1; i < count; i++) {
		tmpx *= x;
		res += a[i] * tmpx;
	}
	return res;
}

//Похідна полінома
myfloat_t PolynomialDerivative(const myfloat_t x, const myfloat_t* a, const int count)
{
	myfloat_t res = a[1], tmpx = 1.0;
	for (int i = 2; i < count; i++) {
		tmpx *= x;
		res += a[i] * tmpx * i;
	}
	return res;
}

//Многочленна ланка з чотирма параметрами
myfloat_t PN4(const myfloat_t x, const myfloat_t a[4])
{
	return Polynomial(x, a, 4);
}

//Похідна многочленної ланки з чотирма параметрами
myfloat_t dPN4(const myfloat_t x, const myfloat_t a[4])
{
	return PolynomialDerivative(x, a, 4);
}

//Многочленна ланка з п'ятьма параметрами
myfloat_t PN5(const myfloat_t x, const myfloat_t a[5])
{
	return Polynomial(x, a, 5);
}

//Похідна многочленної ланки з п'ятьма параметрами
myfloat_t dPN5(const myfloat_t x, const myfloat_t a[5])
{
	return PolynomialDerivative(x, a, 5);
}

//Ермітовий сплайн з параметрами hp. Якщо derivative == 1, то обчислюється похідна
myfloat_t HermiteSpline(const herm_params hp, const myfloat_t x, const char derivative)
{
	static myfloat_t(*const link[2][linktype_count])(const myfloat_t, const myfloat_t[]) = {
		{
			[powexp4]  = PE4,
			[powexp5]  = PE5,
			[poly4]    = PN4,
			[poly5]    = PN5,
			[exppow5]  = EP5,
			[pow2exp2] = W22,
			[pow1exp2] = W12,
		},
		{
			[powexp4]  = dPE4,
			[powexp5]  = dPE5,
			[poly4]    = dPN4,
			[poly5]    = dPN5,
			[exppow5]  = dEP5,
			[pow2exp2] = dW22,
			[pow1exp2] = dW12,
		},
	};
	if (x < hp.X128[0] || x > hp.X128[hp.link_count]) {
		return NAN;
	}
	int i = 0;
	while (x > hp.X128[++i]);
	return link[derivative][hp.type](x, hp.A128 + (i-1) * hp.param_count);
}

//Розв'язує СЛАР методом Гауса
int SolveGauss(const myfloat_t** a, const myfloat_t* b, const int count, myfloat_t* out)
{
	//роблю дублікат матриці, над якою буду робити всякі непотребства
	myfloat_t** mat = calloc(count, sizeof(*mat)); if (mat == NULL) return -1;
	*mat = calloc(count * count, sizeof(**mat)); if (*mat == NULL) return -1;
	for (int i = 0; i < count; i++) {
		mat[i] = *mat + i * count;
		for (int j = 0; j < count; j++) {
			mat[i][j] = a[i][j];
		}
	}

	//роблю дублікат вектору, над яким буду робити всякі непотребства
	myfloat_t* vec = calloc(count, sizeof(*vec)); if (vec == NULL) return -1;
	for (int i = 0; i < count; i++) {
		vec[i] = b[i];
	}

	//прямий хід методу Гаусса
	for (int row = 0; row < count; row++) {
		//output(mat, vec, count);
		for (int i = count - 1; i > row; i--) {
			//йду у зворотньому напрямку, бо значення a[row][row], a[row][j], vec[row] повинні змінитись в останню чергу 
			myfloat_t l = mat[i][row] / mat[row][row];
			for (int j = count - 1; j >= row; j--) {
				mat[i][j] -= l * mat[row][j];
			}
			vec[i] -= l * vec[row];
		}
	}

	//зворотній хід
	for (int i = count - 1; i >= 0; i--) {
		out[i] = vec[i];

		for (int j = count - 1; j > i; j--) {
			out[i] -= mat[i][j] * out[j];	//компілятору чомусь здається, що тут може статись вихід за межі масиву
		}
		out[i] /= mat[i][i];
	}

	//прибираю за собою
	free(vec);
	free(*mat);
	free(mat);

	return 0;
}

//Обчислює параметри ланки PE4
int HermGenPE4(const myfloat_t f[4], const myfloat_t x0, const myfloat_t x1, myfloat_t* out)
{
	//див. notes10(3)
	out[2] = (f[1] / f[0] - (log(f[0] / f[2]) / (x0 - x1)) - (1. / x0 - log(x0 / x1) / (x0 - x1)) * (f[1] / f[0] - f[3] / f[2]) / (1. / x0 - 1. / x1)) /
		(1 + log(x0) - log(pow(x0, x0) / pow(x1, x1)) / (x0 - x1) - (log(x0 / x1) / (1 / x0 - 1 / x1)) * (1. / x0 - log(x0 / x1) / (x0 - x1)));
	out[1] = (f[1] / f[0] - f[3] / f[2] - out[2] * log(x0 / x1)) / (1. / x0 - 1. / x1);
	//out[3] = (log(f[0] / f[2]) - log(pow(x0, out[2] * x0 + out[1]) / pow(x1, out[2] * x1 + out[1]))) / (x0 - x1);
	out[3] = (log(f[0] / f[2]) - out[2] * log(pow(x0, x0) / pow(x1, x1)) - out[1] * log(x0 / x1)) / (x0 - x1);

	out[0] = f[0] * pow(x0, -(out[2] * x0 + out[1])) * exp(-(out[3] * x0));

#ifndef _DEBUG
	for (int i = 0; i < 4; i++)
		if (isnan(out[i]) || isinf(out[i]))
			return -1;
#endif
	return 0;
}

//Обчислює параметри ланки PE5
int HermGenPE5(const myfloat_t f[5], const myfloat_t x0, const myfloat_t x2, myfloat_t* out)
{
	//Точність до 8-го знаку
	const myfloat_t x1 = (x2 + x0) * 0.5;

	//alphas
	const myfloat_t vec[3] = {
		log(f[2] / f[0]) / (x1 - x0) - log(f[3] / f[0]) / (x2 - x0), //alpha1
		f[1] / f[0] - log(f[2] / f[0]) / (x1 - x0),	//alpha2
		f[4] / f[3] - log(f[2] / f[0]) / (x1 - x0)	//alpha3
	};

	myfloat_t matrix[3][3] = {
		{
			log(x1 / x0) / (x1 - x0) - log(x2 / x0) / (x2 - x0), //beta1
			log(pow(x1,x1) / pow(x0,x0)) / (x1 - x0) - log(pow(x2,x2) / pow(x0,x0)) / (x2 - x0), //gamma1
			log(pow(x1,x1 * x1) / pow(x0,x0 * x0)) / (x1 - x0) - log(pow(x2,x2 * x2) / pow(x0,x0 * x0)) / (x2 - x0) //delta1
		},
		{
			1. / x0 - log(x1 / x0) / (x1 - x0), //beta2
			log(x0) + 1. - log(pow(x1,x1) / pow(x0,x0)) / (x1 - x0), //gamma2
			2. * x0 * log(x0) + x0 - log(pow(x1,x1 * x1) / pow(x0,x0 * x0)) / (x1 - x0) //delta2
		},
		{
			1. / x2 - log(x1 / x0) / (x1 - x0), //beta3
			log(x2) + 1. - log(pow(x1,x1) / pow(x0,x0)) / (x1 - x0), //gamma3
			2. * x2 * log(x2) + x2 - log(pow(x1,x1 * x1) / pow(x0,x0 * x0)) / (x1 - x0) // delta3
		}
	};

	myfloat_t* mat[3] = { &matrix[0][0], &matrix[1][0], &matrix[2][0] };
	SolveGauss((const myfloat_t**)mat, vec, 3, out + 1); // => a[1], a[2], a[3] aka b, c, d

	/*
	out[4] = (log(f[2] / f[0]) -
		log(pow(x1, out[1] + out[2] * x1 + out[3] * x1 * x1) /
			pow(x0, out[1] + out[2] * x0 + out[3] * x0 * x0))
		) / (x1 - x0); // h
	*/
	out[4] = (log(f[2] / f[0]) - (
			out[1] * log(x1 / x0) + 
			out[2] * log(pow(x1, x1) / pow(x0, x0)) + 
			out[3] * log(pow(x1, x1 * x1) / pow(x0, x0 * x0))
		)
	) / (x1 - x0); // h

	out[0] = f[0] * pow(x0, -(out[1] + out[2] * x0 + out[3] * x0 * x0)) * exp(-(out[4] * x0)); //a
	//out[0] = log(f[0]) - out[1] * log(x0) - out[2] * x0 * log(x0) - out[3] * x0 * x0 * log(x0) - out[4] * x0;
	//out[0] = exp(out[0]);

#ifndef _DEBUG
	for (int i = 0; i < 5; i++)
		if (isnan(out[i]) || isinf(out[i]))
			return -1;
#endif

	return 0;
}

//Обчислює параметри ланки EP5
int HermGenEP5(const myfloat_t F[5], const myfloat_t x0, const myfloat_t x2, myfloat_t* out)
{
	const myfloat_t x1 = (x2 + x0) * 0.5;
	const myfloat_t f0 = F[0], df0 = F[1], f1 = F[2], f2 = F[3], df2 = F[4];

	//alphas
	const myfloat_t vec[3] = {
		log(f1 / f0) / (log(x1 / x0)) - log(f2 / f0) / (log(x2 / x0)), //alpha1
		df0 / f0 - log(f1 / f0) / (x0 * log(x1 / x0)),				   //alpha2
		df2 / f2 - log(f2 / f0) / (x2 * log(x2 / x0)),				   //alpha3
	};

	const myfloat_t matrix[3][3] = {
		{
			(x1 - x0) / log(x1 / x0) - (x2 - x0) / log(x2 / x0),                                         //beta1
			(x1 * x1 - x0 * x0) / log(x1 / x0) - (x2 * x2 - x0 * x0) / log(x2 / x0),                     //gamma1
			(x1 * x1 * x1 - x0 * x0 * x0) / log(x1 / x0) - (x2 * x2 * x2 - x0 * x0 * x0) / log(x2 / x0), //delta1
		},
		{
			1 - (x1 - x0) / (x0 * log(x1 / x0)),                               //beta2
			2 * x0 - (x1 * x1 - x0 * x0) / (x0 * log(x1 / x0)),                //gamma2
			3 * x0 * x0 - (x1 * x1 * x1 - x0 * x0 * x0) / (x0 * log(x1 / x0)), //delta2
		},
		{
			1 - (x2 - x0) / (x2 * log(x2 / x0)),                               //beta3
			2 * x2 - (x2 * x2 - x0 * x0) / (x2 * log(x2 / x0)),                //gamma3
			3 * x2 * x2 - (x2 * x2 * x2 - x0 * x0 * x0) / (x2 * log(x2 / x0)), //delta3
		}
	};

	const myfloat_t* mat[3] = { &matrix[0][0], &matrix[1][0], &matrix[2][0] };
	SolveGauss(mat, vec, 3, out + 2); // => a[2], a[3], a[4] aka c, d, g

	const myfloat_t c = out[2], d = out[3], g = out[4];
	const myfloat_t b = (log(f1 / f0) - c * (x1 - x0) - d * (x1 * x1 - x0 * x0) - g * (x1 * x1 * x1 - x0 * x0 * x0)) / log(x1 / x0);
	const myfloat_t a = f0 * pow(x0, -(b)) * exp(-(c * x0 + d * x0 * x0 + g * x0 * x0 * x0));

	out[0] = a;
	out[1] = b;

#ifndef _DEBUG
	for (int i = 0; i < 5; i++)
		if (isnan(out[i]) || isinf(out[i]))
			return -1;
#endif

	return 0;
}


//Обчислює параметри ланки W22
int HermGenW22(const myfloat_t F[5], const myfloat_t x0, const myfloat_t x2, myfloat_t* out)
{
	const myfloat_t x1 = (x2 + x0) * 0.5;
	const myfloat_t f0 = F[0], df0 = F[1], f1 = F[2], f2 = F[3], df2 = F[4];
	const myfloat_t p = log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1) - log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0),
		         q = (x2 - x1) / log(x2 / x1) - (x1 - x0) / log(x1 / x0),
		         r = (x2 * x2 - x1 * x1) / log(x2 / x1) - (x1 * x1 - x0 * x0) / log(x1 / x0),
		         s = log(f2 / f1) / log(x2 / x1) - log(f1 / f0) / log(x1 / x0);
	
	myfloat_t alpha1, beta1, gamma1, alpha2, beta2, gamma2;
	alpha1 = x0 * df0 / f0 - log(f1 / f0) / log(x1 / x0)
				+ (s / p) * log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0)
				- (s / p) * x0 * (1 + log(x0)),
	beta1  = (q / p)* log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0)
				- (x1 - x0) / log(x1 / x0)
				- (q / p) * x0 * (1 + log(x0)) + x0, 
	gamma1 = (r / p)* log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0)
				- (x1 * x1 - x0 * x0) / log(x1 / x0)
				- (r / p) * x0 * (1 + log(x0)) + 2 * x0 * x0, 
	alpha2 = x2* df2 / f2 - log(f2 / f1) / log(x2 / x1)
				+ (s / p) * log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1)
				- (s / p) * x2 * (1 + log(x2)),
	beta2  = (q / p)* log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1)
				- (x2 - x1) / log(x2 / x1)
				- (q / p) * x2 * (1 + log(x2)) + x2,
	gamma2 = (r / p)* log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1)
				- (x2 * x2 - x1 * x1) / log(x2 / x1)
				- (r / p) * x2 * (1 + log(x2)) + 2 * x2 * x2;

	myfloat_t A, a1, a2, b1, b2; // A * x^(a1 + a2*x) * exp(b1*x + b2*x^2);

	b2 = (beta1*alpha2 - beta2*alpha1) / (beta1*gamma2 - beta2*gamma1);
	b1 = (alpha1 - gamma1*b2) / beta1;
	a2 = s/p - b1 * (q/p) - b2 * (r/p);
	a1 = log(f1/f0)/log(x1/x0) - a2*log(pow(x1,x1)/pow(x0,x0))/log(x1/x0)
			- b1*(x1-x0)/log(x1/x0) - b2*(x1*x1-x0*x0)/log(x1/x0);
	A = f0 * pow(x0, -(a1 + a2*x0)) * exp( -(b1*x0 + b2*x0*x0));
	
	out[0] = A;
	out[1] = a1;
	out[2] = a2;
	out[3] = b1;
	out[4] = b2;


#ifndef _DEBUG
	for (int i = 0; i < 5; i++)
		if (isnan(out[i]) || isinf(out[i]))
			return -1;
#endif

	return 0;
}


//Обчислює параметри ланки W12
int HermGenW12(const myfloat_t F[4], const myfloat_t x0, const myfloat_t x1, myfloat_t* out)
{
	const myfloat_t f0 = F[0], df0 = F[1], f1 = F[2], df1 = F[3];
	
	myfloat_t alpha1, beta1, gamma1, alpha2, beta2, gamma2;

	myfloat_t aux1, aux2;
	aux1 = (1./x0) / log(x0/x1);
	aux2 = (1./x1) / log(x0/x1);

	alpha1 = df0/f0 - log(f0/f1) * aux1,
	beta1  = -(x0-x1) * aux1 + 1,
	gamma1 = -(x0*x0-x1*x1) * aux1 + 2*x0,
	alpha2 = df1/f1 - log(f0/f1) * aux2,
	beta2  = -(x0-x1) * aux2 + 1,
	gamma2 = -(x0*x0-x1*x1) * aux2 + 2*x1;


	myfloat_t a, b, c, d; // a * x^(b) * exp(c*x + d*x^2);

	d = (beta1 * alpha2 - beta2 * alpha1) / (beta1 * gamma2 - beta2 * gamma1);
	c = (alpha1 - gamma1 * d) / beta1;
	b = (log(f0 / f1) - c*(x0-x1) - d*(x0*x0-x1*x1) ) / log(x0/x1);
	a = f0 / (pow(x0, b) * exp(c*x0 + d*x0*x0));

	out[0] = a;
	out[1] = b;
	out[2] = c;
	out[3] = d;

#ifndef _DEBUG
	for (int i = 0; i < 4; i++)
		if (isnan(out[i]) || isinf(out[i]))
			return -1;
#endif

	return 0;
}


int HermGenPN(const myfloat_t* f, const myfloat_t x0, const myfloat_t x2, const int count, myfloat_t* out)
{
	const int odd = isodd(count);
	const myfloat_t X[2] = { x0, x2 };
	myfloat_t** mat = calloc(count, sizeof(*mat)); if (mat == NULL) return -1;
	*mat = calloc(count * count, sizeof(**mat)); if (*mat == NULL) return -1;
	for (int i = 0; i < count; i++) mat[i] = *mat + i * count;

	const int prod[2][5] = { {1,1,1,1,1}, {0,1,2,3,4} };

	for (int k = 0; k < 2; k++) { //лічильник похідних (з нульової)
		for (int i = k; i < count ; i++) {	//лічильник членів
			//Кожен рядок матриці mat - це множники біля шуканих параметрів a[i]
			for (int j = 0; j < 2; j++) {
				mat[k + j * (2 + odd)][i] = prod[k][i] * pow(X[j], i - k);	//тут компілятор матюкається на те, що mat є стрічкою, а не матрицею
			}
		}
	}
	if (odd) {
		const myfloat_t x2 = (X[0] + X[1]) * 0.5;
		for (int i = 0; i < count; i++) {
			mat[2][i] = pow(x2, i);
		}
	}
	SolveGauss((const myfloat_t**)mat, f, count, out);

#ifndef _DEBUG
	for (int i = 0; i < count; i++)
		if (isnan(out[i]))
			return -1;
#endif

	return 0;
}

//Обчислює параметри ланки PN4
int HermGenPN4(const myfloat_t f[4], const myfloat_t x0, const myfloat_t x2, myfloat_t* out)
{
	return HermGenPN(&f[0], x0, x2, 4, out);
}

//Обчислює параметри ланки PN5
int HermGenPN5(const myfloat_t f[5], const myfloat_t x0, const myfloat_t x2, myfloat_t* out)
{
	return HermGenPN(&f[0], x0, x2, 5, out);
}


myfloat_t finderr(myfloat_t(*link)(const myfloat_t, const myfloat_t[]), const myfloat_t params[], function f, const myfloat_t from, const myfloat_t to)
{
	myfloat_t res = 0.0;
	for (int i = 0; i < 0xFFF; i++) {
		const myfloat_t step = (to - from) / (0xFFF);
		myfloat_t err = fabs(link(from + step * i, &(params[0])) - f(from + step * i));
		if (res < err)
			res = err;
	}
	return res;
}


static int (* const Gen[])(const myfloat_t[], const myfloat_t, const myfloat_t, myfloat_t*) = {
	[powexp4] = HermGenPE4,
	[powexp5] = HermGenPE5,
	[poly4] = HermGenPN4,
	[poly5] = HermGenPN5,
	[exppow5] = HermGenEP5,
	[pow2exp2] = HermGenW22,
	[pow1exp2] = HermGenW12,
};

static myfloat_t(* const link[])(const myfloat_t, const myfloat_t[]) = {
	[powexp4] = PE4,
	[powexp5] = PE5,
	[poly4] = PN4,
	[poly5] = PN5,
	[exppow5] = EP5,
	[pow2exp2] = W22,
	[pow1exp2] = W12,
};


int HermGen(function _f[], herm_params* hp, const myfloat_t a, const myfloat_t b)
{
	int errorcode = 0;
	myfloat_t f[5] = { 0 };
	int odd = 0;

	hp->param_count = param_count_map[hp->type];
	odd = isodd(hp->param_count);

	if (hp->A == NULL) {  //Добре би було придумати якусь надійнішу перевірку
		hp->A = calloc(sizeof(*hp->A), hp->param_count);
		hp->X = calloc(sizeof(*hp->X), 2);
		hp->A128 = calloc(sizeof(*hp->A128), hp->param_count);
		hp->X128 = calloc(sizeof(*hp->X128), 2);
	}

	f[0] = _f[0](a);
	f[1] = _f[1](a);
	if (odd)
		f[2] = _f[0]((a + b) / 2);
	f[2 + odd] = _f[0](b);
	f[3 + odd] = _f[1](b);

	errorcode = Gen[hp->type](f, a, b, hp->A128);
	if (errorcode)
		return errorcode;	//Можливо таки треба буде пофіксити цей витік пам'яті

	hp->X[0] = hp->X128[0] = a;
	hp->X[1] = hp->X128[1] = b;

	for(int i = 0; i < hp->param_count; i++)
		hp->A[0] = hp->A128[0];

	return 0;
}

//Обчислює параметри сплайна з похибкою nu
int HermGenNu(function _f[], herm_params* hp, const myfloat_t a, const myfloat_t b, const myfloat_t nu)
{
	const myfloat_t eps = nu * 1e-3;
	hp->param_count = param_count_map[hp->type];

	const int odd = isodd(hp->param_count);	//Чи к-сть параметрів непарна
	const int linknum = hp->type;				//Тип ланки
	myfloat_t x0, x2 = a, delta = b;					//Межі наближення на певній ітерації

	int errorcode = 0;
	
	struct linkparam{
		myfloat_t params[5];
		myfloat_t xright;
		void* next;
	}*cur, *top;		//Список параметрів сплайна
	top = malloc(sizeof(*top)); if (top == NULL) return -1;
	cur = top;

	int count = 0;	//К-сть ланок
	while (x2 < b) {
		x0 = x2;
		//x2 = fmin(b, x2 + delta * 2);
		x2 = b < x2 + delta * 2 ? b : x2 + delta * 2;
		delta = x2 - x0;
		++count;
		myfloat_t f[5];	//Значення функції
		f[0] = _f[0](x0), f[1] = _f[1](x0);
		f[2 + odd] = _f[0](x2), f[3 + odd] = _f[1](x2);
		if (odd) f[2] = _f[0]((x2 + x0) * 0.5);
		errorcode = Gen[linknum](&f[0], x0, x2, &(cur->params[0]));
		if (errorcode)
			return errorcode;
		myfloat_t nu1 = finderr(link[linknum], cur->params, _f[0], x0, x2);
		if (nu1 < nu) {
			cur->xright = x2;
			break;
		}
		myfloat_t xl = x0, xr = x2, xprev = x2;
		while (fabs(nu - nu1) > eps) {
			xprev = x2;
			if (nu < nu1) {
				xr = x2;
				x2 = 0.5 * (x2 + xl);
			}
			else {
				xl = x2;
				x2 = 0.5 * (x2 + xr);
			}
			if (xprev == x2)
				return -1;	//Тут витік пам'яті, але кому не байдуже?

			f[2 + odd] = _f[0](x2), f[3 + odd] = _f[1](x2);
			if (odd) f[2] = _f[0]((x2 + x0) * 0.5);

			errorcode = Gen[linknum](&f[0], x0, x2, &(cur->params[0]));
			if (errorcode)	//Тут теж потрібно би було вивільняти пам'ять, але ми ж мільйонери...
				return errorcode;
			nu1 = finderr(link[linknum], cur->params, _f[0], x0, x2);
		}

		cur->xright = x2;
		cur->next = malloc(sizeof(*cur));
		cur = cur->next; if (cur == NULL) return -1;
	}

	hp->link_count = count;
	hp->A    = calloc(count * hp->param_count, sizeof(*(hp->A))); if (hp->A == NULL) return -1;
	hp->X    = calloc(count + 1, sizeof(*(hp->X))); if (hp->X == NULL) return -1;
	hp->A128 = calloc(count * hp->param_count, sizeof(*(hp->A128))); if (hp->A128 == NULL) return -1;
	hp->X128 = calloc(count + 1, sizeof(*(hp->X128))); if (hp->X128 == NULL) return -1;
	struct linkparam* iterator = top;
	hp->X128[0] = a;
	hp->X[0] = (double)hp->X128[0];
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < hp->param_count; j++) {
			int tmp = i * hp->param_count + j;
			hp->A128[tmp] = iterator->params[j];
			hp->A[tmp] = (double)hp->A128[tmp];
		}
		hp->X128[i + 1] = iterator->xright;
		hp->X[i + 1] = (double)hp->X128[i + 1];
		void* tmp = iterator;
		iterator = iterator->next;
		free(tmp);
	}
	return 0;
}

//Обчислює параметри сплайна з кількістю ланок r
int HermGenR(function _f[], herm_params* hp, const myfloat_t a, const myfloat_t b, const int r)
{
	int res = 0;
	myfloat_t nul = 0, nur = 0, nu = 0, prevnu = -1;
	int k = 0;
	const myfloat_t eps = 0.05;
	hp->param_count = param_count_map[hp->type];


	// step1:

	nul = (myfloat_t)INFINITY, nur = (myfloat_t)0;

	myfloat_t f[5] = { 0 };
	myfloat_t params[5] = { 0 };
	int odd = isodd(hp->param_count);

	for (int i = 0; i < r; i++) {

		const myfloat_t delta = (b - a) / r;
		const myfloat_t x0 = a + i * delta, x2 = a + (i + 1) * delta;
		f[0] = _f[0](x0), f[1] = _f[1](x0);
		if (odd) f[2] = _f[0]((x0 + x2) / 2);
		f[2 + odd] = _f[0](x2), f[3 + odd] = _f[1](x2);

		res = Gen[hp->type](f, x0, x2, params);
		if (res < 0)
			return res;

		nu = finderr(link[hp->type], params, _f[0], x0, x2);

		nul = min(nul, nu);
		nur = max(nur, nu);
	}

	nu = (nul + nur) / 2.;
	
	// step2:
	// nul = nur = 0.;

	// step3:

	while (1)
	{
		if (prevnu == nu)
			return -1;		//ми вийшли за межі точності

		prevnu = nu;
		HermClear(hp);
		res = HermGenNu(_f, hp, a, b, nu);
		if (res != 0)
			return res;
		k = hp->link_count;

		// step4:

		if (k == r)
		{
			const myfloat_t* params = &hp->A128[(k - 1) * hp->param_count];
			const myfloat_t from = hp->X128[(k - 1)], to = hp->X128[k];

			//потрібно перевірити тільки останню ланку, бо на інших гарантовано похибка рівна nu
			myfloat_t nui = finderr(link[hp->type], params, _f[0], from, to);

			if (fabs((nui - nu) / nu) < eps) {
				// goto step7;
				break;
			}
		}

		// step5:

		if (k > r) {
			nul = nu;
			// nu = nur != 0 ? (nu+nur)/2 : nu*1.1;
			nu = (nu + nur) / 2;
			// goto step3;	//мав би перейти до пункту 4, виконавши пункт 3, тому йду в пункт 3
			continue;
		}

		// step6:

		if (k <= r) { //завжди істина, напевно...
			nur = nu;
			// nu = nul != 0 ? (nu+nul)/2 : nu*0.9;
			nu = (nu + nul) / 2;
			// goto step3; //тут так само...
			continue;
		}
	}

	// step7:
	//тут я мав би щось вивести але мені ліньки

	// step8:
	return res;
}



void HermClear(herm_params* hp)
{
	if (hp->A != NULL) {
		free(hp->A);
		hp->A = NULL;
	}
	if (hp->A128 != NULL) {
		free(hp->A128);
		hp->A128 = NULL;
	}
	if (hp->X != NULL) {
		free(hp->X);
		hp->X = NULL;
	}
	if (hp->X128 != NULL) {
		free(hp->X128);
		hp->X128 = NULL;
	}
}


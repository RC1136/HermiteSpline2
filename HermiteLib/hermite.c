#include "hermite.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <stdio.h>


#define isodd(x) (x & 0b1)

//Степенево-експоненціальна ланка з чотирма параметрами
double PE4(const double x, const double a[4])
{
	return a[0] * pow(x, a[1] + a[2] * x) * exp(a[3] * x);
}

//Похідна степенево-експоненціальної ланки з чотирма параметрами
double dPE4(const double x, const double a[4])
{
	return PE4(x, a) * (a[2] * (1 + log(x)) + a[1] / x + a[3]);
}

//Степенево-експоненціальна ланка з п'ятьма параметрами
double PE5(const double x, const double a[5])
{
	return a[0] * pow(x, a[1] + a[2] * x + a[3] * x * x) * exp(a[4] * x);
}

//Похідна степенево-експоненціальної ланки з п'ятьма параметрами
double dPE5(const double x, const double a[5])
{
	return PE5(x, a) * (
		(a[2] + 2 * a[3] * x) * log(x) + a[1] / x + a[2] + a[3] * x + a[4]
	);
}


double EP5(const double x, const double A[5])
{
	const double a = A[0], b = A[1], c = A[2], d = A[3], g = A[4];
	return a * pow(x, b) * exp(c * x + d * x * x + g * x * x * x);
}


double dEP5(const double x, const double A[5])
{
	const double a = A[0], b = A[1], c = A[2], d = A[3], g = A[4];
	return EP5(x, A) * (b / x + c + 2 * d * x + 3 * g * x * x);
}


double W22(const double x, const double A[5])
{
	// A = A[0], a1 = A[1], a2 = A[2], b1 = A[3], b2 = A[4]
	return A[0] * pow(x, A[1] + A[2] * x) * exp(A[3] * x + A[4] * x * x);
}


double dW22(const double x, const double A[5])
{
	// A = A[0], a1 = A[1], a2 = A[2], b1 = A[3], b2 = A[4]
	return A[0] * pow(x, A[1] + A[2] * x - 1.) * exp(A[3] * x + A[4] * x * x) *
		(A[1] + A[2] * x * (1 + log(x)) + A[3] * x + 2. * A[4] * x * x);
}


//Поліном з count+1 параметрами 
double Polynomial(const double x, const double* a, const int count) 
{
	double res = a[0], tmpx = 1.0;
	for (int i = 1; i < count; i++) {
		tmpx *= x;
		res += a[i] * tmpx;
	}
	return res;
}

//Похідна полінома
double PolynomialDerivative(const double x, const double* a, const int count)
{
	double res = a[1], tmpx = 1.0;
	for (int i = 2; i < count; i++) {
		tmpx *= x;
		res += a[i] * tmpx * i;
	}
	return res;
}

//Многочленна ланка з чотирма параметрами
double PN4(const double x, const double a[4])
{
	return Polynomial(x, a, 4);
}

//Похідна многочленної ланки з чотирма параметрами
double dPN4(const double x, const double a[4])
{
	return PolynomialDerivative(x, a, 4);
}

//Многочленна ланка з п'ятьма параметрами
double PN5(const double x, const double a[5])
{
	return Polynomial(x, a, 5);
}

//Похідна многочленної ланки з п'ятьма параметрами
double dPN5(const double x, const double a[5])
{
	return PolynomialDerivative(x, a, 5);
}

//Ермітовий сплайн з параметрами hp. Якщо derivative == 1, то обчислюється похідна
double HermiteSpline(const herm_params hp, const double x, const char derivative)
{
	static double (*const link[2][linktype_count])(const double, const double[]) = {
		{
			[powexp4]  = PE4,
			[powexp5]  = PE5,
			[poly4]    = PN4,
			[poly5]    = PN5,
			[exppow5]  = EP5,
			[pow2exp2] = W22,
		},
		{
			[powexp4]  = dPE4,
			[powexp5]  = dPE5,
			[poly4]    = dPN4,
			[poly5]    = dPN5,
			[exppow5]  = dEP5,
			[pow2exp2] = dW22,
		},
	};
	if (x < hp.X[0] || x > hp.X[hp.link_count]) {
		return NAN;
	}
	int i = 0;
	while (x > hp.X[++i]);
	return link[derivative][hp.type](x, hp.A + (i-1) * hp.param_count);
}

//Розв'язує СЛАР методом Гауса
int SolveGauss(const double** a, const double* b, const int count, double* out)
{
	//роблю дублікат матриці, над якою буду робити всякі непотребства
	double** mat = calloc(count, sizeof(*mat)); if (mat == NULL) return -1;
	*mat = calloc(count * count, sizeof(**mat)); if (*mat == NULL) return -1;
	for (int i = 0; i < count; i++) {
		mat[i] = *mat + i * count;
		for (int j = 0; j < count; j++) {
			mat[i][j] = a[i][j];
		}
	}

	//роблю дублікат вектору, над яким буду робити всякі непотребства
	double* vec = calloc(count, sizeof(*vec)); if (vec == NULL) return -1;
	for (int i = 0; i < count; i++) {
		vec[i] = b[i];
	}

	//прямий хід методу Гаусса
	for (int row = 0; row < count; row++) {
		//output(mat, vec, count);
		for (int i = count - 1; i > row; i--) {
			//йду у зворотньому напрямку, бо значення a[row][row], a[row][j], vec[row] повинні змінитись в останню чергу 
			double l = mat[i][row] / mat[row][row];
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
int HermGenPE4(const double f[4], const double x0, const double x1, double* out)
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
int HermGenPE5(const double f[5], const double x0, const double x2, double* out)
{
	//Точність до 8-го знаку
	const double x1 = (x2 + x0) * 0.5;

	//alphas
	const double vec[3] = {
		log(f[2] / f[0]) / (x1 - x0) - log(f[3] / f[0]) / (x2 - x0), //alpha1
		f[1] / f[0] - log(f[2] / f[0]) / (x1 - x0),	//alpha2
		f[4] / f[3] - log(f[2] / f[0]) / (x1 - x0)	//alpha3
	};

	double matrix[3][3] = {
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

	double* mat[3] = { &matrix[0][0], &matrix[1][0], &matrix[2][0] };
	SolveGauss(mat, vec, 3, out + 1); // => a[1], a[2], a[3] aka b, c, d

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
int HermGenEP5(const double F[5], const double x0, const double x2, double* out)
{
	const double x1 = (x2 + x0) * 0.5;
	const double f0 = F[0], df0 = F[1], f1 = F[2], f2 = F[3], df2 = F[4];

	//alphas
	const double vec[3] = {
		log(f1 / f0) / (log(x1 / x0)) - log(f2 / f0) / (log(x2 / x0)), //alpha1
		df0 / f0 - log(f1 / f0) / (x0 * log(x1 / x0)),				   //alpha2
		df2 / f2 - log(f2 / f0) / (x2 * log(x2 / x0)),				   //alpha3
	};

	const double matrix[3][3] = {
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

	const double* mat[3] = { &matrix[0][0], &matrix[1][0], &matrix[2][0] };
	SolveGauss(mat, vec, 3, out + 2); // => a[2], a[3], a[4] aka c, d, g

	const double c = out[2], d = out[3], g = out[4];
	const double b = (log(f1 / f0) - c * (x1 - x0) - d * (x1 * x1 - x0 * x0) - g * (x1 * x1 * x1 - x0 * x0 * x0)) / log(x1 / x0);
	const double a = f0 * pow(x0, -(b)) * exp(-(c * x0 + d * x0 * x0 + g * x0 * x0 * x0));

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
int HermGenW22(const double F[5], const double x0, const double x2, double* out)
{
	const double x1 = (x2 + x0) * 0.5;
	const double f0 = F[0], df0 = F[1], f1 = F[2], f2 = F[3], df2 = F[4];
	const double p = log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1) - log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0),
		         q = (x2 - x1) / log(x2 / x1) - (x1 - x0) / log(x1 / x0),
		         r = (x2 * x2 - x1 * x1) / log(x2 / x1) - (x1 * x1 - x0 * x0) / log(x1 / x0),
		         s = log(f2 / f1) / log(x2 / x1) - log(f1 / f0) / log(x1 / x0);
	
	double alpha1, beta1, gamma1, alpha2, beta2, gamma2;
	alpha1 = x0 * df0 / f0 - log(f1 / f0) / f(x1 / x0)
				+ (s / p) * log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0)
				- (s / p) * x0 * (1 + log(x0)),
	beta1  = (q / p)* log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0)
				- (x1 - x0) / log(x1 / x0)
				- (q / p) * x0 * (1 + log(x0)) + x0, 
	gamma1 = (r / p)* log(pow(x1, x1) / pow(x0, x0)) / log(x1 / x0)
				- (x1 * x1 - x0 * x0) / log(x1 / x0)
				- (r / p) * x0 * (1 + log(x0)) + 2 * x0 * x0, 
	alpha2 = x2* df2 / f2 - log(f2 / f1) / f(x2 / x1)
				+ (s / p) * log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1)
				- (s / p) * x2 * (1 + log(x2)),
	beta2  = (q / p)* log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1)
				- (x2 - x1) / log(x2 / x1)
				- (q / p) * x2 * (1 + log(x2)) + x2,
	gamma2 = (r / p)* log(pow(x2, x2) / pow(x1, x1)) / log(x2 / x1)
				- (x2 * x2 - x1 * x1) / log(x2 / x1)
				- (r / p) * x2 * (1 + log(x2)) + 2 * x2 * x2;

	double A, a1, a2, b1, b2; // A * x^(a1 + a2*x) * exp(b1*x + b2*x^2);

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


int HermGenPN(const double* f, const double x0, const double x2, const int count, double* out)
{
	const int odd = isodd(count);
	const double X[2] = { x0, x2 };
	double** mat = calloc(count, sizeof(*mat)); if (mat == NULL) return -1;
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
		const double x2 = (X[0] + X[1]) * 0.5;
		for (int i = 0; i < count; i++) {
			mat[2][i] = pow(x2, i);
		}
	}
	SolveGauss(mat, f, count, out);

	for (int i = 0; i < 4; i++)
		if (isnan(out[i]))
			return -1;

	return 0;
}

//Обчислює параметри ланки PN4
int HermGenPN4(const double f[4], const double x0, const double x2, double* out)
{
	return HermGenPN(&f[0], x0, x2, 4, out);
}

//Обчислює параметри ланки PN5
int HermGenPN5(const double f[5], const double x0, const double x2, double* out)
{
	return HermGenPN(&f[0], x0, x2, 5, out);
}


double finderr(double (*link)(const double, const double[]), const double params[], function f, const double from, const double to)
{
	double res = 0.0;
	for (int i = 0; i < 0xFFF; i++) {
		const double step = (to - from) / (0xFFF);
		double err = fabs(link(from + step * i, &(params[0])) - f(from + step * i));
		if (res < err)
			res = err;
	}
	return res;
}


//Обчислює параметри сплайна з похибкою nu
int HermGen(function _f[], herm_params* hp, const double a, const double b, const double nu)
{
	// FILE* logger = fopen("log.txt", "w+");
	// fprintf(logger, "1\r\n");
	// fflush(logger);
	const double eps = nu * 1e-5;
	static int (* const Gen[])(const double[], const double, const double, double*) = {
		[powexp4]  = HermGenPE4,
		[powexp5]  = HermGenPE5,
		[poly4]    = HermGenPN4,
		[poly5]    = HermGenPN5,
		[exppow5]  = HermGenEP5,
		[pow2exp2] = HermGenW22,
	};
	static double (* const link[])(const double, const double[]) = { 
		[powexp4]  = PE4,
		[powexp5]  = PE5,
		[poly4]    = PN4,
		[exppow5]  = EP5,
		[poly5]    = PN5,
		[pow2exp2] = W22,
	};
	// fprintf(logger, "2\r\n");
	// fflush(logger);

	const int odd = isodd(hp->param_count),	//Чи к-сть параметрів непарна
		linknum = hp->type;				//Тип ланки
	double x0, x2 = a, delta = b;					//Межі наближення на певній ітерації

	int errorcode = 0;
	

	// fprintf(logger, "3\r\n");
	// fflush(logger);
	struct linkparam{
		double params[5];
		double xright;
		void* next;
	}*cur, *top;		//Список параметрів сплайна
	top = malloc(sizeof(*top)); if (top == NULL) return -1;
	cur = top;

	// fprintf(logger, "4\r\n");
	// fflush(logger);
	int count = 0;	//К-сть ланок
	while (x2 < b) {
		// fprintf(logger, "5\r\n");
		// fflush(logger);
		x0 = x2;
		x2 = fmin(b, x2 + delta * 2);
		delta = x2 - x0;
		++count;
		double f[5];	//Значення функції
		f[0] = _f[0](x0), f[1] = _f[1](x0);
		f[2 + odd] = _f[0](x2), f[3 + odd] = _f[1](x2);
		if (odd) f[2] = _f[0]((x2 + x0) * 0.5);
		errorcode = Gen[linknum](&f[0], x0, x2, &(cur->params[0]));
		if (errorcode)
			return errorcode;
		double nu1 = finderr(link[linknum], cur->params, _f[0], x0, x2);
		if (nu1 < nu) {
			cur->xright = x2;
			break;
		}
		// fprintf(logger, "6\r\n");
		// fflush(logger);
		double xl = x0, xr = x2, xprev = x2;
		while (fabs(nu - nu1) > eps) {
			// fprintf(logger, "% .30lf\t% .30lf\t% .30lf\t% .30lf\t% .30le\n", x0, xl, x2, xr, nu1);
			// fflush(logger);
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
			
			for (int i = 0; i < hp->param_count; i++)
				printf("%lf\n", f[i]);
			putchar('\n');
			
			errorcode = Gen[linknum](&f[0], x0, x2, &(cur->params[0]));
			if (errorcode)	//Тут теж потрібно би було вивільняти пам'ять, але ми ж мільйонери...
				return errorcode;
			nu1 = finderr(link[linknum], cur->params, _f[0], x0, x2);
		}
		/*
		for (int i = 0; i < hp->param_count; i++) {
			printf("A[%d][%d] == %lf\n", count - 1, i, cur->params[i]);
		}
		putchar('\n');
		putchar('\n');
		*/
		cur->xright = x2;
		cur->next = malloc(sizeof(*cur));
		cur = cur->next; if (cur == NULL) return -1;
	}

	hp->link_count = count;
	hp->A = calloc(count * hp->param_count, sizeof(*(hp->A))); if (hp->A == NULL) return -1;
	hp->X = calloc(count + 1, sizeof(*(hp->X))); if (hp->X == NULL) return -1;
	struct linkparam* iterator = top;
	hp->X[0] = a;
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < hp->param_count; j++) {
			hp->A[i * hp->param_count + j] = iterator->params[j];
			//printf("A[%d][%d] == %lf\n", i, j, iterator->params[j]);
		}
		hp->X[i + 1] = iterator->xright;
		void* tmp = iterator;
		iterator = iterator->next;
		free(tmp);
	}
	//putchar('\n');

	// fclose(logger);
	return 0;
}

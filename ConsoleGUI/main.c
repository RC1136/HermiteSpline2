
#include <stdio.h>
#include <math.h>

#include <Windows.h> // дл€ SetConsoleOutputCP

#include "HermiteLib/hermite.h"
#include "HermiteLib/functions.h"

#define countof(x) ( sizeof(x) / sizeof(*(x)) )

int main(void)
{
	SetConsoleOutputCP(1251);

	int linknum = 0;
	printf("1. a[0] * x^(a[1] + a[2]*x) * exp(a[3]*x)              \n");
	printf("2. a[0] * x^(a[1] + a[2]*x + a[3]*x^2) * exp(a[4]*x)   \n");
	printf("3. a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3                 \n");
	printf("4. a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3 + a[4]*x^4      \n");
	printf("5. a[0] + x^(a[1]) * exp(a[2]*x + a[3]*x^2 + a[4]*x^3) \n");

	printf("¬ибер≥ть вигл€д ланки сплайна (введ≥ть число в≥д 1 до 5): ");

	int res = scanf("%d", &linknum);
	if (res < 1)
		return -1;

	int funcnum = 0;
	printf("0.  1x^(4x - 3) * exp(7x)         \n");
	printf("1.  1x^(-2 + 5x + 3x^2) * exp(2x) \n");
	printf("2.  3 + 2x - 4x^2 + x^3           \n");
	printf("3.  3 + 2x - 4x^2 + x^3 - x^4     \n");
	printf("4.  exp( (x - 4)^3 / 42 )         \n");
	printf("5.  sin(x)                        \n");
	printf("6.  sqrt(x^3 + 1)                 \n");
	printf("7.  1 / (x^2 + 0.4)               \n");
	printf("8.  tg(x) / (x + 2) + 2           \n");
	printf("9.  ln(x)                         \n");
	printf("10. exp( sin(x) - cos(x) )        \n");
	printf("11. exp(sin(x)) + x               \n");
	printf("12. tg(x)                         \n");
	printf("13. (-3 * x^2) / (x^3 + 1)^2      \n");

	printf("¬ибер≥ть наближувану функц≥ю (введ≥ть число в≥д 0 до 13): ");

	res = scanf("%d", &funcnum);
	if (res < 1)
		return -1;

	double a = 0., b = 0.;
	printf("¬вед≥ть меж≥ наближенн€: ");
	res = scanf("%lf%lf", &a, &b);
	if (res < 2)
		return -1;

	double nu = 0.;
	printf("¬вед≥ть точн≥сть наближенн€: ");
	res = scanf("%lf", &nu);

	if (res == 0)
		nu = INFINITY;

	herm_params hp = { 0 };
	hp.type = linknum;

	// ѕерепиши це, будь ласка
	if ((linknum == 1) || (linknum == 3))
		hp.param_count = 4;
	else
		hp.param_count = 5;

	function f[] = { funcs[funcnum], dfuncs[funcnum]};
	HermGen(f, &hp, a, b, nu);


	return 0;
}
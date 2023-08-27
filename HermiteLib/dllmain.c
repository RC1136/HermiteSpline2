#include <windows.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "hermite.h"
#include "functions.h"

herm_params __declspec(dllexport) _HermGen(int8_t funcnum, int8_t linknum, double a, double b, double nu)
{
	herm_params hp = { 0 };
	hp.type = linknum;

	function f[] = { funcs[funcnum], dfuncs[funcnum] };
	HermGen(f, &hp, a, b, nu);

	return hp;
}

void __declspec(dllexport) _free(herm_params hp)
{
	free(hp.A);
	free(hp.X);
}

double __declspec(dllexport) _HermiteSpline(const herm_params hp, const double x, int8_t der)
{
	return HermiteSpline(hp, x, der);
}

double __declspec(dllexport) _Func(const int8_t funcnum, const double x, int8_t der)
{
	return der ? dfuncs[funcnum](x) : funcs[funcnum](x);
}

double __declspec(dllexport) _MaxError(const herm_params hp, const int8_t funcnum, const double from, const double to)
{
	double res = 0.0;
	for (int i = 0; i < 0xFFFF; i++) {
		const double step = (to - from) / (0xFFFF);
		double err = fabs(HermiteSpline(hp, from + step * i, 0) - funcs[funcnum](from + step * i));
		if (res < err)
			res = err;
	}
	return res;
}

BOOL APIENTRY DllMain(HMODULE hModule,
	DWORD ul_reason_for_call,
	LPVOID lpReserved) {
	switch (ul_reason_for_call) {
	case DLL_PROCESS_ATTACH:
		break;
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}










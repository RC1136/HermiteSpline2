#include <windows.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "float_def.h"
#include "hermite.h"
#include "functions.h"

herm_params __declspec(dllexport) _HermGenNu(int8_t funcnum, int8_t linknum, double a, double b, double nu)
{
	herm_params hp = { 0 };
	hp.type = linknum;

	function f[] = { funcs[funcnum], dfuncs[funcnum] };
	int err = HermGenNu(f, &hp, (myfloat_t)a, (myfloat_t)b, (myfloat_t)nu);
	if (err != 0) {
		hp.type = 0;
	}
	return hp;
}

herm_params __declspec(dllexport) _HermGenR(int8_t funcnum, int8_t linknum, double a, double b, int32_t linkcount)
{
	herm_params hp = { 0 };
	hp.type = linknum;

	function f[] = { funcs[funcnum], dfuncs[funcnum] };
	int err = HermGenR(f, &hp, (myfloat_t)a, (myfloat_t)b, (int)linkcount);
	if (err != 0) {
		hp.type = 0;
	}
	return hp;
}

void __declspec(dllexport) _free(herm_params hp)
{
	HermClear(&hp);
}

double __declspec(dllexport) _HermiteSpline(const herm_params hp, const double x, int8_t der)
{
	if (hp.type == 0)
	{
		return -42.;
	}
	return (double)HermiteSpline(hp, (myfloat_t)x, der);
}

double __declspec(dllexport) _Func(const int8_t funcnum, const double x, int8_t der)
{
	return (double)(der ? dfuncs[funcnum]((myfloat_t)x) : funcs[funcnum]((myfloat_t)x));
}

double __declspec(dllexport) _MaxError(const herm_params hp, const int8_t funcnum, const double from, const double to)
{
	if (hp.type == 0)
	{
		return -42.;
	}
	myfloat_t res = 0.0;
	const myfloat_t step = ((myfloat_t)to - (myfloat_t)from) / ((myfloat_t)0xFFFF);
	for (int i = 0; i < 0xFFFF; i++) {
		myfloat_t err = fabs(HermiteSpline(hp, (myfloat_t)from + step * i, 0) - funcs[funcnum]((myfloat_t)from + step * i));
		if (res < err)
			res = err;
	}
	return (double)res;
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










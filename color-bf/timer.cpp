/*#include <Windows.h>

long long milliseconds_now()
{
	static LARGE_INTEGER s_frequency;
	static bool s_use_qpc = QueryPerformanceFrequency(&s_frequency);
	if (s_use_qpc) {
		LARGE_INTEGER now;
		QueryPerformanceCounter(&now);
		return (1000LL * now.QuadPart) / s_frequency.QuadPart;
	} else {
		return GetTickCount();
	}
}*/
#include <iostream>
using namespace std;

long long milliseconds_now()
{
    return clock();
}
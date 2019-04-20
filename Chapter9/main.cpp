#include "IncompressCouetteFlow.h"
#include <fstream>

#ifdef _DEBUG
#pragma comment(lib, "CFDMath_d.lib")
#else
#pragma comment(lib, "CFDMath.lib")
#endif

int main()
{
	double E1 = 1;
	double ReD = 5000;
	size_t N = 20;
	size_t IterateTimes1 = 360;
	std::set<size_t> SelectedIterateTimes1 = { 12, 36, 60, 120, 240, 360 };
	std::ofstream fout1("ImplicitIncompressCouetteFlow.csv");
	ImplicitIncompressCouetteFlow(E1, ReD, N, IterateTimes1, SelectedIterateTimes1, fout1);
	double E2 = 0.1;
	size_t IterateTimes2 = 1000;
	std::set<size_t> SelectedIterateTimes2 = { 12, 36, 60, 120, 240, 360, 600, 800, 1000 };
	std::ofstream fout2("ExplicitIncompressCouetteFlow.csv");
	ExplicitIncompressCouetteFlow(E2, ReD, N, IterateTimes2, SelectedIterateTimes2, fout2);
	return 0;
}
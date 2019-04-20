#include "../CFDMath/CFDMath.h"
#include "IncompressCouetteFlow.h"

using namespace std;
using namespace CFDMath;

void ImplicitIncompressCouetteFlow(double E,
	double ReD,
	size_t N,
	size_t IterateTimes,
	set<size_t> SelectedIterateTimes,
	ostream &out)
{
	assert(N > 2);
	double delta_y = 1.0 / N;
	double delta_t = E * ReD * delta_y * delta_y;
	double A = -E / 2;
	double B = 1 + E;

	vector<double> u(N - 1, 0.0), K(N - 1, 0.0);
	vector<double> a(N - 2, A), b(N - 1, B), c(N - 2, A);
	for (size_t n = 1; n <= IterateTimes; ++n)
	{
		K[0] = (1 - E) * u[0] + E * u[1] / 2;
		for (size_t j = 1; j < N - 2; ++j)
		{
			K[j] = (1 - E) * u[j] + E * (u[j - 1] + u[j + 1]) / 2;
		}
		K[N - 2] = (1 - E) * u[N - 2] + E * (u[N - 3] + 1) / 2 - A;
		TridiagonalMatrixThomas(a, b, c, K, u);
		if (SelectedIterateTimes.find(n) != SelectedIterateTimes.end())
		{
			out << n << "," << 0.0 << ",";
			for (size_t j = 0; j < N - 1; ++j)
			{
				out << u[j] << ",";
			}
			out << 1.0 <<"\n";
		}
	}
}
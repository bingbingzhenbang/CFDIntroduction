#include "../CFDMath/CFDMath.h"
#include "IncompressCouetteFlow.h"

using namespace std;
using namespace CFDMath;

void ExplicitIncompressCouetteFlow(double E,
	double ReD,
	size_t N,
	size_t IterateTimes,
	set<size_t> SelectedIterateTimes,
	ostream &out)
{
	assert(N > 2);
	double delta_y = 1.0 / N;
	double delta_t = E * ReD * delta_y * delta_y;
	vector<double> u(N - 1, 0.0), v(N - 1, 0.0);
	for (size_t n = 1; n <= IterateTimes; ++n)
	{
		v[0] = u[0] + E * (u[1] - 2 * u[0]);
		for (size_t j = 1; j < N - 2; ++j)
		{
			v[j] = u[j] + E * (u[j + 1] - 2 * u[j] + u[j - 1]);
		}
		v[N - 2] = u[N - 2] + E * (1 - 2 * u[N - 2] + u[N - 3]);
		u.swap(v);
		if (SelectedIterateTimes.find(n) != SelectedIterateTimes.end())
		{
			out << n << "," << 0.0 << ",";
			for (size_t j = 0; j < N - 1; ++j)
			{
				out << u[j] << ",";
			}
			out << 1.0 << "\n";
		}
	}
}
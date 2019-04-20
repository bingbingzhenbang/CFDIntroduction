#include "CFDMath.h"

using namespace std;

namespace CFDMath {

	void TridiagonalMatrixThomas(const vector<double> &a,
		const vector<double> &b,
		const vector<double> &c,
		const vector<double> &d,
		vector<double> &x)
	{
		if (a.empty() || b.empty() || c.empty() || d.empty())
			return;
		assert(a.size() == c.size());
		assert(b.size() == d.size());
		assert(a.size() + 1 == b.size());

		size_t N = b.size();
		x.swap(vector<double>(N, 0.0));

		vector<double> l(N - 1, 0.0);
		vector<double> u(N, 0.0);
		u[0] = b[0];
		for (size_t i = 1; i < N; ++i)
		{
			l[i - 1] = a[i - 1] / u[i - 1];
			u[i] = b[i] - l[i - 1] * c[i - 1];
		}

		vector<double> y(N, 0.0);
		y[0] = d[0];
		for (size_t i = 1; i < N; ++i)
		{
			y[i] = d[i] - l[i - 1] * y[i - 1];
		}
		x[N - 1] = y[N - 1] / u[N - 1];
		for (size_t i = N - 1; i > 0; --i)
		{
			x[i - 1] = (y[i - 1] - c[i - 1] * x[i]) / u[i - 1];
		}
	}

}
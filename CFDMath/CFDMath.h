#ifndef CFDMATH_H
#define CFDMATH_H

#include <cassert>
#include <vector>

namespace CFDMath {

	void TridiagonalMatrixThomas(const std::vector<double> &a,
		const std::vector<double> &b,
		const std::vector<double> &c,
		const std::vector<double> &d,
		std::vector<double> &x);

}

#endif // CFDMATH_H
#ifndef ABSR_BSPLINE_H
#define ABSR_BSPLINE_H

#include <common.h>

namespace absr {
	class BSplines {
	public:
		BSplines() {};

		BSplines(const Vector &_knots);

		Scalar B(const int &i, const int &k, const Scalar &t) const;

		Scalar dB(const int &i, const int &k, const Scalar &t) const;

		Scalar d2B(const int &i, const int &k, const Scalar &t) const;

		Scalar w(const int &i, const int &k, const Scalar &t) const;

		Scalar X(const int &i, const Scalar &t) const;

	protected:

		Vector knots_;
	};
};

#endif
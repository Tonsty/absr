#include <bspline.h>

using namespace absr;

BSplines::BSplines(const Vector &_knots) : knots_(_knots) {}

Scalar BSplines::B(const int &i, const int &k, const Scalar &t) const {
	if(t >= knots_(i+k) || t < knots_(i)) return 0.0;
	if(k == 1) return X(i, t);
	return w(i, k, t) * B(i, k-1, t) + (1 - w(i+1, k, t)) * B(i+1, k-1, t);
}

Scalar BSplines::dB(const int &i, const int &k, const Scalar &t) const {
	if(k <= 1) return 0.0;

	Scalar rt1 = 0.0, rt2 = 0.0;
	if(knots_(i+k) != knots_(i+1)) rt1 = -B(i+1, k-1, t) / (knots_(i+k) - knots_(i+1));
	if(knots_(i+k-1) != knots_(i)) rt2 = B(i, k-1, t) / (knots_(i+k-1) - knots_(i));

	return (k-1) * (rt1 + rt2);
}

Scalar BSplines::d2B(const int &i, const int &k, const Scalar &t) const {
	if(k <= 2) return 0.0;

	Scalar rt1 = 0.0, rt2 = 0.0;
	if(knots_(i+k) != knots_(i+1)) rt1 = -dB(i+1, k-1, t) / (knots_(i+k) - knots_(i+1));
	if(knots_(i+k-1) != knots_(i)) rt2 = dB(i, k-1, t) / (knots_(i+k-1) - knots_(i));

	return (k-1) * ( rt1 + rt2);
}

Scalar BSplines::w(const int &i, const int &k, const Scalar &t) const {
	if (knots_(i) != knots_(i+k-1)) {
		return (t - knots_(i)) / (knots_(i+k-1) - knots_(i));
	}
	return 0.0;
}

Scalar BSplines::X(const int &i, const Scalar &t) const {
	if(t >= knots_(i) && t < knots_(i+1)) return 1;
	return 0.0;
}
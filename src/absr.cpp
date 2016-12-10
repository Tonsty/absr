#include <absr.h>

using namespace absr;

void ABSR::build_sdf() {}

void ABSR::abspline_fitting_sdf(Scalar lambda) {
	tensorbs_.fitting_sdf(sdf_, lambda);
}

void ABSR::abspline_fitting_3L(Scalar lambda) {
	tensorbs_.fitting_3L(points_, normals_, lambda);
}

void ABSR::abspline_fitting_Juttler(Scalar lambda, Scalar kappa) {
	tensorbs_.fitting_Juttler(points_, normals_, lambda, kappa);
}

void ABSR::resample_sdf(SDF &sdf) {
	PointSet points;
	sdf.topoints(points);
	tensorbs_.evaluate(points, sdf.values_);
}
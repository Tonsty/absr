#include <absr.h>

using namespace absr;

void ABSR::build_sdf() {}

void ABSR::abspline_fitting_sdf() {
	tensorbs_.fitting_sdf(sdf_);
}

void ABSR::abspline_fitting_3L() {
	tensorbs_.fitting_3L(points_, normals_);
}

void ABSR::resample_sdf(SDF &sdf) {
	PointSet points;
	sdf.topoints(points);
	tensorbs_.evaluate(points, sdf.values_);
}
#ifndef ABSR_TENSORBSPLINES_H
#define ABSR_TENSORBSPLINES_H

#include <common.h>
#include <sdf.h>

namespace absr {
	struct TensorBSplines {
		TensorBSplines() {};
		void fitting_sdf(const SDF &sdf, Scalar lambda = 0.1);
		void fitting_3L(const PointSet &points, const NormalSet &normals, Scalar lambda = 0.07);
		void fitting_Juttler(const PointSet &points, const NormalSet &normals, Scalar lambda = 0.08, Scalar kappa = 0.05);
		void evaluate(const PointSet &points, Vector &values);
		Vector controls_;
	};
};

#endif

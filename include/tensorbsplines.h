#ifndef ABSR_TENSORBSPLINES_H
#define ABSR_TENSORBSPLINES_H

#include <common.h>
#include <sdf.h>

namespace absr {
	struct TensorBSplines {
		TensorBSplines() {};
		void fitting_sdf(const SDF &sdf);
		void fitting_3L(const PointSet &points, const NormalSet &normals);
		void fitting_Juttler(const PointSet &points, const NormalSet &normals);
		void evaluate(const PointSet &points, Vector &values);
		Vector controls_;
	};
};

#endif

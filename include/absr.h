#ifndef ABSR_ABSR_H
#define ABSR_ABSR_H

#include <common.h>
#include <sdf.h>
#include <tensorbsplines.h>

namespace absr {

	struct Mesh {
		PointSet verts_;
		Faces faces_;
	};

	struct BoundingBox {
		Point low_;
		Point high_;
	};

	class ABSR {

	public:

		ABSR() {};

		ABSR(const SDF &sdf) : sdf_(sdf) {}

		ABSR(const PointSet &points, const NormalSet &normals) : points_(points), normals_(normals) {};

		void compute_boundingbox();

		void compute_normals();

		void build_sdf();

		void abspline_fitting_sdf();

		void abspline_fitting_3L();

		void abspline_fitting_Juttler();

		void resample_sdf(SDF &sdf);

	protected:

		TensorBSplines tensorbs_;

		PointSet points_;
		NormalSet normals_;
		SDF sdf_;
	};
}

#endif
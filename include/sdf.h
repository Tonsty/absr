#ifndef ABSR_SDF_H
#define ABSR_SDF_H

#include <typedefs.h>

namespace absr {
	struct SDF {
		SDF() {}
		void topoints_1d(Vector &points_1d) const;
		void topoints(PointSet &points) const;
		Vector values_;
		Size grid_size_;
		Scalar voxel_length_;
	};
}

#endif
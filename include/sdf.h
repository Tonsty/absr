#ifndef ABSR_SDF_H
#define ABSR_SDF_H

#include <common.h>

namespace absr {
	struct SDF {
		SDF() {}
		void topoints(PointSet &points) const;
		Vector values_;
		Size grid_size_;
		Scalar voxel_length_;
	};
}

#endif
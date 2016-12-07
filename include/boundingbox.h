#ifndef ABSR_BOUNDINGBOX_H
#define ABSR_BOUNDINGBOX_H

#include <common.h>

namespace absr {
	struct BoundingBox {
		Point low_;
		Point high_;
		BoundingBox(const PointSet &points);
		static TransformMat normalize_to_unit_cube(PointSet &points);
	};
};

#endif
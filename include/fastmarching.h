#ifndef ABSR_FASTMARCHING_H
#define ABSR_FASTMARCHING_H

#include <sdf.h>

namespace absr {

	struct FastMarching {
		void compute(const PointSet &points, Scalar narrow_band_width = (Scalar) 3.5);
		void tagging();
		Size grid_size_; 
		Scalar voxel_length_;
		Vector values_;
		enum State {INACTIVE, ACTIVE, FIXED};
		std::vector<State> flag_;
	};

};

#endif
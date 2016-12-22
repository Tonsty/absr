#ifndef POLYGONIZER_HPP
#define POLYGONIZER_HPP

#include <tensorbsplines.h>

namespace absr {
	struct Polygonizer {
		void compute(TensorBSplines &tbs, Size grid_size, Scalar isovalue = 0.0,
			TransformMat transmat = TransformMat::Identity(4,4), bool invertface = false);
	};
};

#endif



#ifndef POLYGONIZER_HPP
#define POLYGONIZER_HPP

#include <tensorbsplines.h>

namespace absr {
	struct Polygonizer {
		void compute(Function *f, Size grid_size, Scalar isovalue = 0.0, Scalar x = 0.5, Scalar y = 0.5, Scalar z = 0.5,
			TransformMat transmat = TransformMat::Identity(4,4), bool invertface = false);
	};
};

#endif



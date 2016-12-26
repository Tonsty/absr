#ifndef ABSR_IO_H
#define ABSR_IO_H

#include <typedefs.h>

namespace absr {
	struct IO {
		static void load_points_normals(const std::string file, PointSet &points, NormalSet &normals);
	};
};

#endif
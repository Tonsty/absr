#ifndef ABSR_IO_H
#define ABSR_IO_H

#include <typedefs.h>

namespace absr {
	struct IO {
		static void load_points_normals(const std::string points_file, const std::string normals_file, PointSet &points, NormalSet &normals);
	};
};

#endif
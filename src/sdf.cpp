#include <sdf.h>

using namespace absr;

void SDF::topoints(PointSet &points) const {
	const Size grid_size = grid_size_;
	const Scalar voxel_length = voxel_length_;

	Size npts = grid_size*grid_size*grid_size;
	points.resize(npts, 3);

	for(Index zi = 0; zi < grid_size; zi++) {
		for(Index yi = 0; yi < grid_size; yi++) {
			for (Index xi = 0; xi < grid_size; xi++) {
				Index row_index = xi+yi*grid_size+zi*grid_size*grid_size;
				Scalar xpos = xi*voxel_length, ypos = yi*voxel_length, zpos = zi*voxel_length;
				points.row(row_index) << xpos, ypos, zpos;
			}
		}
	}
}
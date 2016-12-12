#include <fastmarching.h>
#include <vector>
#include <queue>
#include <iostream>

using namespace absr;

//struct QNode {
//	QNode() {}
//	QNode(Index _index, Scalar _distance) : index_(_index), distance_(_distance) {}
//	bool operator<(const QNode&other) const {
//		return distance_ < other.distance_;
//	}
//	Index index_;
//	Scalar distance_;
//};

struct RNode {
	RNode() {}
	RNode(Index _index, Scalar _distance) : index_(_index), distance_(_distance) {}
	bool operator<(const RNode&other) const {
		return distance_ > other.distance_;
	}
	Index index_;
	Scalar distance_;
};


void FastMarching::compute(const PointSet &points, Scalar narrow_band_width) {

	const Size grid_size = grid_size_;
	const Scalar voxel_length = voxel_length_;

	Scalar MAX_DIST = narrow_band_width*voxel_length;//std::numeric_limits<Scalar>::max();

	std::cerr << "\nbegin fast marching:" << std::endl;
	std::cerr << "grid = "<< grid_size << " * " << grid_size << " * " << grid_size  << ", (MAX_DIST = " << MAX_DIST <<" )" << std::endl;

	values_ = Vector::Constant(grid_size*grid_size*grid_size, MAX_DIST);
	flag_ = std::vector<State>(grid_size*grid_size*grid_size, INACTIVE);

	Size npts = points.rows();
	std::vector<RNode> initial_interface;
	for (Index pi = 0; pi < npts; pi++) {
		Point pt = points.row(pi);
		Scalar xpos = pt(0), ypos = pt(1), zpos = pt(2);
		Index xi = (Index)(xpos/voxel_length), yi = (Index)(ypos/voxel_length), zi = (Index)(zpos/voxel_length);

		Index xbegin = ( (xpos-xi*voxel_length)/voxel_length > (Scalar) 0.5 ) ? xi : ((xi-1)>=0?(xi-1):0);
		Index ybegin = ( (ypos-yi*voxel_length)/voxel_length > (Scalar) 0.5 ) ? yi : ((yi-1)>=0?(yi-1):0);
		Index zbegin = ( (zpos-zi*voxel_length)/voxel_length > (Scalar) 0.5 ) ? zi : ((zi-1)>=0?(zi-1):0);

		Index xend = (xbegin+2) < grid_size ? (xbegin+2) : (grid_size-1);
		Index yend = (ybegin+2) < grid_size ? (ybegin+2) : (grid_size-1);
		Index zend = (zbegin+2) < grid_size ? (zbegin+2) : (grid_size-1);

		for (Index zii = zbegin; zii <= zend; zii++) {
			for (Index yii = ybegin; yii <= yend; yii++) {
				for (Index xii = xbegin; xii <= xend; xii++) {
					Index gi = xii + yii * grid_size + zii * grid_size * grid_size;
					if (flag_[gi]==INACTIVE) {
						initial_interface.push_back(RNode(gi, (Scalar)(-1.0)));
						flag_[gi] = ACTIVE; }
					Point grid_pt(3);
					grid_pt << xii*voxel_length, yii*voxel_length, zii*voxel_length;
					values_(gi) = std::min(values_(gi), (pt-grid_pt).norm());
				}
			}
		}
	}

	std::priority_queue<RNode> pq;
	for (Index i = 0; i < initial_interface.size(); i++) {
		initial_interface[i].distance_ = values_(initial_interface[i].index_);
		pq.push(initial_interface[i]);
	}

	while (!pq.empty()) {
		RNode t =pq.top();
		Index index = t.index_;
		Scalar distance = t.distance_;
		pq.pop();

		//if (values_(index) < distance) continue;
		flag_[index] = FIXED;

		Index zi = index/(grid_size*grid_size);
		index -= zi*grid_size*grid_size;
		Index yi = index/grid_size;
		index -= yi*grid_size;
		Index xi = index;

		Point pt(3);
		pt << xi*voxel_length, yi*voxel_length, zi*voxel_length;
		
		Index xbegin = (xi-1) >= 0 ? (xi-1) : 0;
		Index ybegin = (yi-1) >= 0 ? (yi-1) : 0;
		Index zbegin = (zi-1) >= 0 ? (zi-1) : 0;

		Index xend = (xi+1) < grid_size ? (xi+1) : (grid_size-1);
		Index yend = (yi+1) < grid_size ? (yi+1) : (grid_size-1);
		Index zend = (zi+1) < grid_size ? (zi+1) : (grid_size-1);

		for (Index zii = zbegin; zii <= zend; zii++) {
			for (Index yii = ybegin; yii <= yend; yii++) {
				for (Index xii = xbegin; xii <= xend; xii++) {
					Index gi = xii + yii * grid_size + zii * grid_size * grid_size;
					if (flag_[gi] == FIXED) continue;
					Point grid_pt(3);
					grid_pt << xii*voxel_length, yii*voxel_length, zii*voxel_length;
					Scalar old_distance = values_(gi);
					Scalar new_distance = distance + (pt-grid_pt).norm();
					grid_pt << xii*voxel_length, yii*voxel_length, zii*voxel_length;
					if (flag_[gi]==INACTIVE && new_distance < MAX_DIST
						//|| new_distance < old_distance
						) {
						pq.push(RNode(gi, new_distance));
						values_(gi) = new_distance;
						flag_[gi] = ACTIVE; 
					}
				}
			}
		}
	}

	std::cerr << "finished fast marching" << std::endl;
}

void FastMarching::tagging() {

	std::cerr << "\nbegin tagging:" << std::endl;

	const Size grid_size = grid_size_;
	const Scalar voxel_length = voxel_length_;

	values_ = -values_; //invert sign
	flag_ = std::vector<State>(grid_size*grid_size*grid_size, INACTIVE);

	std::vector<RNode> initial_interface;
	for (Index yi = 0; yi < grid_size; yi++) {
		for (Index xi = 0; xi < grid_size; xi++) {
			Index gi = yi*grid_size + xi;
			if (flag_[gi]==INACTIVE) {
				initial_interface.push_back(RNode(gi, values_(gi)));
				flag_[gi] = ACTIVE;
			}
		}
	}
	for (Index zi = 0; zi < grid_size; zi++) {
		for (Index xi = 0; xi < grid_size; xi++) {
			Index gi = zi*grid_size*grid_size + xi;
			if (flag_[gi]==INACTIVE) {
				initial_interface.push_back(RNode(gi, values_(gi)));
				flag_[gi] = ACTIVE;
			}
		}
	}
	for (Index zi = 0; zi < grid_size; zi++) {
		for (Index yi = 0; yi < grid_size; yi++) {
			Index gi = zi*grid_size*grid_size + yi*grid_size;
			if (flag_[gi]==INACTIVE) {
				initial_interface.push_back(RNode(gi, values_(gi)));
				flag_[gi] = ACTIVE;
			}
		}
	}

	for (Index yi = 0; yi < grid_size; yi++) {
		for (Index xi = 0; xi < grid_size; xi++) {
			Index gi = (grid_size-1)*grid_size*grid_size + yi*grid_size+xi;
			if (flag_[gi]==INACTIVE) {
				initial_interface.push_back(RNode(gi, values_(gi)));
				flag_[gi] = ACTIVE;
			}
		}
	}
	for (Index zi = 0; zi < grid_size; zi++) {
		for (Index xi = 0; xi < grid_size; xi++) {
			Index gi = zi*grid_size*grid_size + (grid_size-1)*grid_size + xi;
			if (flag_[gi]==INACTIVE) {
				initial_interface.push_back(RNode(gi, values_(gi)));
				flag_[gi] = ACTIVE;
			}
		}
	}
	for (Index zi = 0; zi < grid_size; zi++) {
		for (Index yi = 0; yi < grid_size; yi++) {
			Index gi = zi*grid_size*grid_size + yi*grid_size + (grid_size-1);
			if (flag_[gi]==INACTIVE) {
				initial_interface.push_back(RNode(gi, values_(gi)));
				flag_[gi] = ACTIVE;
			}
		}
	}

	std::priority_queue<RNode> pq;
	for (Index i = 0; i < initial_interface.size(); i++) pq.push(initial_interface[i]);	

	while (!pq.empty()) {
		RNode t =pq.top();
		Index index = t.index_;
		Scalar distance = t.distance_;
		pq.pop();

		if(abs(distance) < voxel_length) break;

		flag_[index] = FIXED;
		values_(index) = -values_(index);

		Index zi = index/(grid_size*grid_size);
		index -= zi*grid_size*grid_size;
		Index yi = index/grid_size;
		index -= yi*grid_size;
		Index xi = index;

		Index xbegin = (xi-1) >= 0 ? (xi-1) : 0;
		Index ybegin = (yi-1) >= 0 ? (yi-1) : 0;
		Index zbegin = (zi-1) >= 0 ? (zi-1) : 0;

		Index xend = (xi+1) < grid_size ? (xi+1) : (grid_size-1);
		Index yend = (yi+1) < grid_size ? (yi+1) : (grid_size-1);
		Index zend = (zi+1) < grid_size ? (zi+1) : (grid_size-1);

		for (Index zii = zbegin; zii <= zend; zii++) {
			for (Index yii = ybegin; yii <= yend; yii++) {
				for (Index xii = xbegin; xii <= xend; xii++) {
					Index gi = xii + yii * grid_size + zii * grid_size * grid_size;
					if (flag_[gi] == FIXED) continue;
					Scalar its_distance = values_(gi);
					if (flag_[gi]==INACTIVE && its_distance >= distance) {
						flag_[gi] = ACTIVE; 
						pq.push(RNode(gi, its_distance));
					}
				}
			}
		}
	}

	std::cerr << "finished tagging" << std::endl;
}
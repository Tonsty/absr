#include <hierarchicaltbs.h>
#include <tensorbsplines.h>
#include <iostream>

namespace absr {
	template<Degree deg>
	void HierarchicalTBS<deg>::evaluate(const PointSet &points, Vector &values) {
		const Size L = L_;
		evaluate_toplevels(points, L, values);
	}

	template<Degree deg>
	void HierarchicalTBS<deg>::evaluate_toplevels(const PointSet &points, const Size toplevels, Vector &values) {
		Size npts = (Size) points.rows();
		values = Vector::Zero(npts);

		std::cerr << "begin evaluation:" << std::endl;

		for (Index row_index = 0; row_index < npts; row_index++) {
			Point point = points.row(row_index);
			Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();
			values(row_index) = (*this)(xpos, ypos, zpos, toplevels);
		}

		std::cerr << "finished evaluation" << std::endl;
	}

	template<Degree deg> 
	Scalar HierarchicalTBS<deg>::operator()(Scalar x, Scalar y, Scalar z) {
		const Size L = L_;
		return (*this)(x, y, z, L);
	}

	template<Degree deg>
	Scalar HierarchicalTBS<deg>::operator()(Scalar x, Scalar y, Scalar z, const Size toplevels) {

		const Size root_N = root_N_;

		Scalar sum_value = 0;
		Index control_begin_index = 0;
		for (Size level = 0; level < toplevels; level++) {

			const MapType &amp = amps_[level]; 
			const Size current_N = (root_N-deg)*(1<<level) + deg;

			std::vector<std::pair<Index, Scalar>> iws;
			TensorBSplines<deg>::make_tbs_iws(iws, x, y, z, current_N);

			for(auto it = iws.begin(); it!= iws.end(); it++) {
				Index virtual_control_index = it->first;
				Scalar weight = it->second;
				Scalar coefficient;
				auto amp_it = amp.find(virtual_control_index);
				if(amp_it == amp.end()) {		
					coefficient = 0.0; //evaluate at inactive region
				} else {
					Index real_control_index = amp_it->second;
					coefficient = controls_(control_begin_index + real_control_index);
				}
				sum_value += coefficient * weight;
			}

			control_begin_index += (Size) amp.size();
		}

		return sum_value;
	}
};

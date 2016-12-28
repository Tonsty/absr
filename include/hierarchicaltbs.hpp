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

	template<Degree deg>
	void HierarchicalTBS<deg>::precompute_ti0ti1ti2(Index bottomlevel, Vector &ti0, Vector &ti1, Vector &ti2) {
		Size tablesize = (deg+1)*(1<<bottomlevel) + deg;
		ti0 = Vector::Zero(tablesize);
		ti1 = Vector::Zero(tablesize);
		ti2 = Vector::Zero(tablesize);

		Function1D *pfs0 = (Function1D*)(&(TensorBSplines<deg>::pfs0));
		ScaleShiftFunction1D pfs0_h(pfs0, (Scalar) 1.0/(1<<bottomlevel));

		Function1D *pfs1 = (Function1D*)(&(TensorBSplines<deg>::pfs1));
		ScaleShiftFunction1D pfs1_h(pfs1, (Scalar) 1.0/(1<<bottomlevel)); 

		Function1D *pfs2 = (Function1D*)(&(TensorBSplines<deg>::pfs2));
		ScaleShiftFunction1D pfs2_h(pfs2, (Scalar) 1.0/(1<<bottomlevel)); 

		for (Index i = 0; i < tablesize; i++) {

			ScaleShiftFunction1D pfs0_l(pfs0, (Scalar) 1.0, (Scalar) deg-i); 
			ProductFunction1D product_hl0((Function1D*)(&pfs0_h), (Function1D*)(&pfs0_l));
			for (Index h = 0; h < deg+1; h++) {
				ti0(i) += integration((Function1D*)(&product_hl0), deg+deg, (Scalar)(i-deg+h), (Scalar)(i-deg+h+1)); 
			}

			ScaleShiftFunction1D pfs1_l(pfs1, (Scalar) 1.0, (Scalar) deg-i); 
			ProductFunction1D product_hl1((Function1D*)(&pfs1_h), (Function1D*)(&pfs1_l));
			for (Index h = 0; h < deg+1; h++) {
				ti1(i) += integration((Function1D*)(&product_hl1), deg-1+deg-1, (Scalar)(i-deg+h), (Scalar)(i-deg+h+1));
			}

			ScaleShiftFunction1D pfs2_l(pfs2, (Scalar) 1.0, (Scalar) deg-i); 
			ProductFunction1D product_hl2((Function1D*)(&pfs2_h), (Function1D*)(&pfs2_l));
			for (Index h = 0; h < deg+1; h++) {
				ti2(i) += integration((Function1D*)(&product_hl2), deg-2+deg-2,(Scalar)(i-deg+h), (Scalar)(i-deg+h+1));
			}
		}
	}

	template<Degree deg>
	static void HierarchicalTBS<deg>::make_smooth_triplelist(const MapVecType &amps, const Index row_level, const Index column_level,
		const Index row_start_index, const Index column_start_index, std::vector<Triplet> &smooth_triplelist, const Size root_N) {

	}

	template<Degree deg>
	void HierarchicalTBS<deg>::make_smooth_mat_vec(const MapVecType &amps, SparseMatrix &global_smooth_mat, const Vector &controls,
		const Index currentlevel, SparseMatrix &smooth_mat, Vector &smooth_vec, 
		std::vector<Vector> &ti0s, std::vector<Vector> &ti1s, std::vector<Vector> &ti2s, const Size root_N) {

			const Size current_N = (root_N-deg)*(1<<currentlevel) + deg;
			const MapType &amp = amps[currentlevel];

			Vector t0, t1, t2;
			HierarchicalTBS<deg>::precompute_ti0ti1ti2(currentlevel, t0, t1, t2);

			if (currentlevel == 0) {
				ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, current_N);
				smooth_vec = Vector::Zero(amp.size());
				global_smooth_mat = smooth_mat;
			} else {
				const SparseMatrix &A = global_smooth_mat;

				std::vector<Triplet> B_tripleList;
				Size B_rows_cnt = 0;
				for (Index level = 0; level < currentlevel-1; level++) {
					const MapType &amp_pre = amps[level];

					B_rows_cnt += (Size) amp_pre.size();
				}
				SparseMatrix B(B_rows_cnt, amp.size());
				B.setFromTriplets(B_tripleList.begin(), B_tripleList.end());

				SparseMatrix C(amp.size(), amp.size());
				ActiveTBS<deg>::make_smooth_mat(amp, C, current_N);

				Eigen::SparseMatrix<Scalar, Eigen::RowMajor> BC(B_rows_cnt + amp.size(), amp.size());
				BC.topRows(B_rows_cnt) = B; BC.bottomRows(amp.size()) = C;

				Eigen::SparseMatrix<Scalar, Eigen::RowMajor> AB(B_rows_cnt + amp.size(), B_rows_cnt);
				AB.topRows(B_rows_cnt) = A; AB.bottomRows(amp.size()) = B.transpose();

				smooth_mat = BC.transpose() * BC;
				smooth_vec = -BC.transpose() * (AB * controls); 
			}
			ti0s.push_back(t0); ti1s.push_back(t1); ti2s.push_back(t2);
	}
};

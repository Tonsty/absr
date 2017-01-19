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
	void HierarchicalTBS<deg>::make_smooth_triplelist(const MapVecType &amps, const Index row_level, const Index column_level, 
		const Index row_start_index, const Index column_start_index, std::vector<Triplet> &smooth_triplelist, 
		const std::vector<Vector> &ti0s, const std::vector<Vector> &ti1s, const std::vector<Vector> &ti2s, const Size root_N) {
		
		const MapType &column_amp = amps[column_level];
		const MapType &row_amp = amps[row_level];

		const Size column_N = ((root_N-deg)<<column_level) + deg;
		const Size row_N = ((root_N-deg)<<row_level) + deg;

		Scalar delta = (Scalar) 1.0/(column_N-deg);
		Scalar delta3 = delta * delta * delta;

		for (auto column_amp_it = column_amp.begin(); column_amp_it != column_amp.end(); column_amp_it++) {
			Index virtual_column_index = column_amp_it->first;
			Index real_column_index = column_amp_it->second;

			Index column_k = virtual_column_index/(column_N*column_N);
			virtual_column_index -= column_k*column_N*column_N;
			Index column_j = virtual_column_index/column_N;
			virtual_column_index -= column_j*column_N;
			Index column_i = virtual_column_index;

			Index level_diff = column_level - row_level;
			Index row_k = (column_k + (deg<<level_diff) - deg) >> level_diff;
			Index row_j = (column_j + (deg<<level_diff) - deg) >> level_diff;
			Index row_i = (column_i + (deg<<level_diff) - deg) >> level_diff;

			for (Index kk = std::max(row_k-deg, 0); kk <= std::min(row_k+deg, row_N-1); kk++) {
				for (Index jj = std::max(row_j-deg, 0); jj <= std::min(row_j+deg, row_N-1); jj++) {
					for (Index ii = std::max(row_i-deg, 0); ii <= std::min(row_i+deg, row_N-1); ii++) {

						Index t_k_index = (column_k + (deg<<level_diff) - deg) - (kk<<level_diff) + deg;
						Index t_j_index = (column_j + (deg<<level_diff) - deg) - (jj<<level_diff) + deg;
						Index t_i_index = (column_i + (deg<<level_diff) - deg) - (ii<<level_diff) + deg;

						if (t_k_index>=0 && t_j_index>=0 && t_i_index>=0) {
							Index virtual_row_index = kk * row_N * row_N + jj * row_N + ii;
							auto row_amp_it = row_amp.find(virtual_row_index);
							if (row_amp_it != row_amp.end()) {
								Index real_row_index = row_amp_it->second;
								Scalar weight = ti2s[level_diff](t_i_index) * ti0s[level_diff](t_j_index) * ti0s[level_diff](t_k_index) 
									+ ti0s[level_diff](t_i_index) * ti2s[level_diff](t_j_index) * ti0s[level_diff](t_k_index) 
									+ ti0s[level_diff](t_i_index) * ti0s[level_diff](t_j_index) * ti2s[level_diff](t_k_index) 
									+ ti1s[level_diff](t_i_index) * ti1s[level_diff](t_j_index) * ti0s[level_diff](t_k_index) * 2
									+ ti1s[level_diff](t_i_index) * ti0s[level_diff](t_j_index) * ti1s[level_diff](t_k_index) * 2
									+ ti0s[level_diff](t_i_index) * ti1s[level_diff](t_j_index) * ti1s[level_diff](t_k_index) * 2;
								weight *= delta3;
								smooth_triplelist.push_back(Triplet(real_row_index + row_start_index, 
									real_column_index + column_start_index, weight));
							}
						}
					}
				}
			}
		}
	}

	template<Degree deg>
	void HierarchicalTBS<deg>::make_smooth_mat_vec(const MapVecType &amps, std::vector<Triplet> &global_smooth_triplelist, const Vector &controls,
		const Index currentlevel, SparseMatrix &smooth_mat, Vector &smooth_vec, 
		const std::vector<Vector> &ti0s, const std::vector<Vector> &ti1s, const std::vector<Vector> &ti2s, const Size root_N) {

			const MapType &amp = amps[currentlevel];

			if (currentlevel == 0) {
				HierarchicalTBS<deg>::make_smooth_triplelist(amps, 0, 0, 0, 0, global_smooth_triplelist, ti0s, ti1s, ti2s, root_N);
				smooth_mat.resize(amp.size(), amp.size());
				smooth_mat.setFromTriplets(global_smooth_triplelist.begin(), global_smooth_triplelist.end());
				smooth_vec = Vector::Zero(amp.size());
			} else {
				SparseMatrix A(controls.size(), controls.size());
				A.setFromTriplets(global_smooth_triplelist.begin(), global_smooth_triplelist.end());

				std::vector<Triplet> B_tripleList;
				Size B_rows_cnt = 0;
				for (Index level = 0; level < currentlevel; level++) {
					HierarchicalTBS<deg>::make_smooth_triplelist(amps, level, currentlevel, B_rows_cnt, 0, B_tripleList, ti0s, ti1s, ti2s, root_N);
					B_rows_cnt += (Size) amps[level].size();
				}
				SparseMatrix B(B_rows_cnt, amp.size());
				B.setFromTriplets(B_tripleList.begin(), B_tripleList.end());

				std::vector<Triplet> C_tripleList;
				HierarchicalTBS<deg>::make_smooth_triplelist(amps, currentlevel, currentlevel, 0, 0, C_tripleList, ti0s, ti1s, ti2s, root_N);
				SparseMatrix C(amp.size(), amp.size());
				C.setFromTriplets(C_tripleList.begin(), C_tripleList.end());

				Eigen::SparseMatrix<Scalar, Eigen::RowMajor> BC(B_rows_cnt + amp.size(), amp.size());
				BC.topRows(B_rows_cnt) = B; BC.bottomRows(amp.size()) = C;

				Eigen::SparseMatrix<Scalar, Eigen::RowMajor> ABt(B_rows_cnt + amp.size(), B_rows_cnt);
				ABt.topRows(B_rows_cnt) = A; ABt.bottomRows(amp.size()) = B.transpose();

				smooth_mat = BC.transpose() * BC;
				smooth_vec = -BC.transpose() * (ABt * controls); 

				std::vector<Triplet> newB_tripleList, newBt_tripleList;
				for (auto B_it = B_tripleList.begin(); B_it != B_tripleList.end(); B_it++) {
					newB_tripleList.push_back(Triplet(B_it->row(), B_it->col() + B_rows_cnt, B_it->value()));
					newBt_tripleList.push_back(Triplet(B_it->col() + B_rows_cnt, B_it->row(), B_it->value()));
				}

				global_smooth_triplelist.insert(global_smooth_triplelist.end(), newB_tripleList.begin(), newB_tripleList.end());
				global_smooth_triplelist.insert(global_smooth_triplelist.end(), newBt_tripleList.begin(), newBt_tripleList.end());

				std::vector<Triplet> newC_tripleList;
				for (auto C_it = C_tripleList.begin(); C_it != C_tripleList.end(); C_it++) {
					newC_tripleList.push_back(Triplet(C_it->row() + B_rows_cnt, C_it->col() + B_rows_cnt, C_it->value()));
				}
				global_smooth_triplelist.insert(global_smooth_triplelist.end(), newC_tripleList.begin(), newC_tripleList.end());
			}
	}
};

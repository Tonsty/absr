#include <activetbs.h>
#include <tensorbsplines.h>
#include <iostream>

namespace absr {
	template <Degree deg>
	void ActiveTBS<deg>::evaluate(const PointSet &points, Vector &values) {
		Size npts = (Size) points.rows();
		values = Vector::Zero(npts);

		std::cerr << "begin evaluation:" << std::endl;

		for (Index row_index = 0; row_index < npts; row_index++) {
			Point point = points.row(row_index);
			Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();
			values(row_index) = (*this)(xpos, ypos, zpos);
		}

		std::cerr << "finished evaluation" << std::endl;
	}

	template <Degree deg>
	Scalar ActiveTBS<deg>::operator()(Scalar x, Scalar y, Scalar z) {
		const Size N = N_; 
		const MapType &amp = amp_; 

		std::vector<std::pair<Index, Scalar>> iws;
		TensorBSplines<deg>::make_tbs_iws(iws, x, y, z, N);

		Scalar sum_value = 0;
		for(auto it = iws.begin(); it!= iws.end(); it++) {
			Index virtual_control_index = it->first;
			Scalar weight = it->second;
			Scalar coefficient;
			auto amp_it = amp.find(virtual_control_index);
			if(amp_it == amp.end()) {		
				coefficient = 1.0; //evaluate at inactive region
			} else {
				Index real_control_index = amp_it->second;
				coefficient = controls_(real_control_index);
			}
			sum_value += coefficient * weight;
		}

		return sum_value;
	}

	template <Degree deg>
	void ActiveTBS<deg>::generate_active_map(MapType &amp, const PointSet &points, const Size N) {
		const Size npts = (Size) points.rows();
		const Scalar delta = (Scalar) 1.0/(N-deg);

		if(N < deg+1) {
			std::cerr << "for deg = " << deg << " b-splines, number of knots should be at least " << deg+1 << std::endl;
			exit(0);
		}else if (N == deg+1) {
			Size N3 = N * N * N;
			for (Index virtual_index = 0; virtual_index < N3; virtual_index++) {
				Index real_index = virtual_index;
				amp[virtual_index] = real_index;
			}
		}else {
			for (Index row_index = 0; row_index < npts; row_index++) {
				Scalar xpos = points(row_index, 0), ypos = points(row_index, 1), zpos = points(row_index, 2);
				Index i = (Index)(xpos/delta), j = (Index)(ypos/delta), k = (Index)(zpos/delta);
				i = (i>=N-deg) ? (N-deg-1) : i;
				j = (j>=N-deg) ? (N-deg-1) : j;
				k = (k>=N-deg) ? (N-deg-1) : k;

				for(Index t = 0; t <= deg; t++) {
					for(Index s = 0; s <= deg; s++) {
						for(Index r = 0; r <= deg; r++) {
							Index virtual_index = i+r+(j+s)*N+(k+t)*N*N;
							auto amp_it = amp.find(virtual_index);
							if(amp_it == amp.end()) {
								Size real_index = (Size) amp.size();
								amp[virtual_index] = real_index;
							}
						}
					}
				}
			}
		}
	}

	template <Degree deg>
	void ActiveTBS<deg>::make_data_mat(const MapType &amp, const PointSet &points, SparseMatrix &data_mat, const Size N) {

		std::cerr << "\nprepare data_mat:" << std::endl;

		const Size npts = (Size) points.rows();
		const Scalar delta = (Scalar) 1.0/(N-deg);

		data_mat.resize(npts, amp.size());

		std::vector<Triplet> data_list;
		for (Index row_index = 0; row_index < npts; row_index++) {
			Point point = points.row(row_index);
			Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();
			IndexWeightVec iws;

			TensorBSplines<deg>::make_tbs_iws(iws, xpos, ypos, zpos, N);
			for(auto it = iws.begin(); it!= iws.end(); it++) {
				Index virtual_column_index = it->first;
				auto amp_it = amp.find(virtual_column_index);
				if(amp_it == amp.end()) {
					std::cerr << "make_data_mat: visit inactive regions" << std::endl;
					exit(0);
				} else {
					Index real_column_index = amp_it->second;
					Scalar weight = it->second;
					data_list.push_back(Triplet(row_index, real_column_index, weight));
				}
			}
		}
		data_mat.setFromTriplets(data_list.begin(), data_list.end());

		std::cerr << "finished data_mat" << std::endl;
	}

	template <Degree deg>
	void ActiveTBS<deg>::make_smooth_mat(const MapType &amp, SparseMatrix &smooth_mat, const Size N) {

		std::cerr << "\nprepare smooth_mat:" << std::endl;

		Matrix Pi0, Pi1, Pi2;
		TensorBSplines<deg>::precompute_Pi0Pi1Pi2(Pi0, Pi1, Pi2, N);

		//std::cerr << "Pi0 = \n" << Pi0 << std::endl;
		//std::cerr << "Pi1 = \n" << Pi1 << std::endl;
		//std::cerr << "Pi2 = \n" << Pi2 << std::endl;

		std::vector<Triplet> smooth_tripleList;
		Size cnt = 0;
		for (auto amp_it = amp.begin(); amp_it != amp.end(); amp_it++) {
			Index real_row_index = amp_it->second;
			Index virtual_row_index = amp_it->first;

			Index k = virtual_row_index/(N*N);
			virtual_row_index -= k*N*N;
			Index j = virtual_row_index/N;
			virtual_row_index -= j*N;
			Index i = virtual_row_index;

			for (Index kk = std::max(k-deg, 0); kk <= std::min(k+deg, N-1); kk++) {
				for (Index jj = std::max(j-deg, 0); jj <= std::min(j+deg, N-1); jj++) {
					for (Index ii = std::max(i-deg, 0); ii <= std::min(i+deg, N-1); ii++) {
						Index virtual_column_index = kk*N*N + jj*N + ii;
						auto amp_it2 = amp.find(virtual_column_index);
						if(amp_it2 != amp.end()) {
							Index real_column_index = amp_it2->second;
							if(real_column_index < real_row_index) continue;
							Scalar weight = Pi2(i, ii) * Pi0(j, jj) * Pi0(k, kk) 
								+ Pi0(i, ii) * Pi2(j, jj) * Pi0(k, kk)
								+ Pi0(i, ii) * Pi0(j, jj) * Pi2(k, kk)
								+ Pi1(i, ii) * Pi1(j, jj) * Pi0(k, kk) * 2
								+ Pi1(i, ii) * Pi0(j, jj) * Pi1(k, kk) * 2
								+ Pi0(i, ii) * Pi1(j, jj) * Pi1(k, kk) * 2;
							smooth_tripleList.push_back(Triplet(real_row_index, real_column_index, weight));
							if(real_column_index > real_row_index) smooth_tripleList.push_back(Triplet(real_column_index, real_row_index, weight));
						}
					}
				}
			}
		}

		smooth_mat.resize(amp.size(), amp.size());
		smooth_mat.setFromTriplets(smooth_tripleList.begin(), smooth_tripleList.end());

		std::cerr << "finished smooth_mat" << std::endl;
	}

	template <Degree deg>
	void ActiveTBS<deg>::make_data_duvw_mat(const MapType &amp, const PointSet &points, SparseMatrix &data_mat, SparseMatrix &du_mat, 
		SparseMatrix &dv_mat, SparseMatrix &dw_mat, const Size N) {

			std::cerr << "prepare data_mat, du_mat, dv_mat, dw_mat" << std::endl;

			Size npts = (Size) points.rows();
			data_mat.resize(npts, amp.size());
			du_mat.resize(npts, amp.size());
			dv_mat.resize(npts, amp.size());
			dw_mat.resize(npts, amp.size());

			std::vector<Triplet> data_list, du_list, dv_list, dw_list;
			for (Index row_index = 0; row_index < npts; row_index++) {
				Point point = points.row(row_index);
				Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();

				IndexWeightVec iws;
				std::vector<Scalar> du_ws, dv_ws, dw_ws;
				TensorBSplines<deg>::make_tbs_iws_duvw_ws(iws, du_ws, dv_ws, dw_ws, xpos, ypos, zpos, N);
				Index h = 0;
				for(auto it = iws.begin(); it!= iws.end(); it++, h++) {
					Index virtual_column_index = it->first;
					auto amp_it = amp.find(virtual_column_index);
					if(amp_it == amp.end()) {
						std::cerr << "make_data_mat: visit inactive regions" << std::endl;
						exit(0);
					} else {
						Index real_column_index = amp_it->second;
						Scalar weight = it->second, du_weight = du_ws[h], dv_weight = dv_ws[h], dw_weight = dw_ws[h];
						data_list.push_back(Triplet(row_index, real_column_index, weight));
						du_list.push_back(Triplet(row_index, real_column_index, du_weight));
						dv_list.push_back(Triplet(row_index, real_column_index, dv_weight));
						dw_list.push_back(Triplet(row_index, real_column_index, dw_weight));
					}
				}
			}
			data_mat.setFromTriplets(data_list.begin(), data_list.end());
			du_mat.setFromTriplets(du_list.begin(), du_list.end());
			dv_mat.setFromTriplets(dv_list.begin(), dv_list.end());
			dw_mat.setFromTriplets(dw_list.begin(), dw_list.end());

			std::cerr << "finished data_mat, du_mat, dv_mat, dw_mat" << std::endl;
	}
};
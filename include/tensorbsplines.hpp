#include <tensorbsplines.h>
#include <iostream>

namespace absr {

	template<Degree deg>
	PolyFormSpline<deg, 0> TensorBSplines<deg>::pfs0;
	template<Degree deg>
	PolyFormSpline<deg, 1> TensorBSplines<deg>::pfs1;
	template<Degree deg>
	PolyFormSpline<deg, 2> TensorBSplines<deg>::pfs2;

	Vector MULT(const Vector &U, const Vector &V);
	Matrix INTofMULT(const Matrix &A, const Matrix &B);

	template<Degree deg, Degree der>
	PolyFormSpline<deg, der>::PolyFormSpline() {
		switch (deg) {
		case 2: {
			B.resize(3, 3);
			switch(der) {
			case 0: { B << 1, -2, 1, 1, 2, -2, 0, 0, 1; break; }
			case 1: { B << -2, 2, 0, 2, -4, 0, 0, 2, 0; break; }
			case 2: { B << 2, 0, 0, -4, 0, 0, 2, 0, 0; break;}
			default: { std::cerr << "der must be 0, 1 or 2" << std::endl; exit(0); } 
			};
			B/=2;
			break;
			};
		case 3: {
			B.resize(4, 4);
			switch(der) {
			case 0: { B<< 1, -3, 3, -1, 4, 0, -6, 3, 1, 3, 3, -3,	0, 0, 0, 1; break;}
			case 1: { B<< -3, 6, -3, 0, 0, -12, 9, 0, 3, 6, -9, 0, 0, 0, 3, 0; break;}
			case 2: { B<< 6, -6, 0, 0, -12, 18, 0, 0, 6, -18, 0, 0, 0, 6, 0, 0; break;}
			default: { std::cerr << "der must be 0, 1 or 2" << std::endl; exit(0);}
			};
			B/=6;	
			break;
			};
		default: { std::cerr << "deg must be 2 or 3" << std::endl; exit(0);}
		}
	}

	template<Degree deg, Degree der>
	Scalar PolyFormSpline<deg, der>::operator()(Scalar x) {
		if(x < 0 || x > (deg+1)) return 0;
		Index i = (Index) x;
		i = i < (deg+1) ? i : deg;
		Scalar u = x - i;
		Vector pu = Vector::Ones(deg+1);
		for (Index h = 0; h < deg; h++) pu(h+1) = pu(h) * u;
		return B.row(deg-i) * pu;
	}

	template<Degree deg>
	Scalar TensorBSplines<deg>::operator()(Scalar x, Scalar y, Scalar z) {
		const Size N = N_; 

		std::vector<std::pair<Index, Scalar>> iws;
		make_tbs_iws(iws, x, y, z, N);

		Scalar sum_value = 0;
		for(auto it = iws.begin(); it!= iws.end(); it++) {
			Index control_index = it->first;
			Scalar weight = it->second;
			sum_value += controls_(control_index) * weight;
		}

		return sum_value;
	}

	template<Degree deg> 
	void TensorBSplines<deg>::evaluate(const PointSet &points, Vector &values) {
		const Size N = N_; 
		Size npts = (Size) points.rows();
		values = Vector::Zero(npts);

		std::cerr << "begin evaluation:" << std::endl;

		Scalar delta = (Scalar) 1.0/(N-deg);
		Arrayi ijk = (points/delta).array().floor().cast<Index>();
		for(Index h = 0; h < npts; h++) ijk(h,0) = ijk(h,0) < (N-deg) ? ijk(h,0) : (N-deg-1);
		for(Index h = 0; h < npts; h++) ijk(h,1) = ijk(h,1) < (N-deg) ? ijk(h,1) : (N-deg-1);
		for(Index h = 0; h < npts; h++) ijk(h,2) = ijk(h,2) < (N-deg) ? ijk(h,2) : (N-deg-1);

		Array uvw = (points/delta).array() - ijk.cast<Scalar>();

		const Matrix &B = pfs0.B;
		Array U = Array::Ones(npts, deg+1), V = Array::Ones(npts, deg+1), W = Array::Ones(npts, deg+1);
		for (Index h = 0; h < deg; h++) {
			U.col(h+1) = U.col(h) * uvw.col(0);
			V.col(h+1) = V.col(h) * uvw.col(1);
			W.col(h+1) = W.col(h) * uvw.col(2);
		}

		Matrix bu = U.matrix()*B.transpose(), bv = V.matrix()*B.transpose(), bw = W.matrix()*B.transpose();

		for (Index kk = 0; kk <= deg; kk++) {
			for (Index jj = 0; jj <= deg; jj++) {
				for (Index ii = 0; ii <= deg; ii++) {
					Arrayi indices = N*N*(ijk.col(2)+kk) + N*(ijk.col(1)+jj) + ijk.col(0)+ii;
					Array bubvbw = bu.col(ii).array()*bv.col(jj).array()*bw.col(kk).array();
					for(Index h = 0; h < npts; h++) values(h) += controls_(indices(h))*bubvbw(h);
				}
			}
		}

		//for (Index row_index = 0; row_index < npts; row_index++) {
		//	Point point = points.row(row_index);
		//	Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();
		//	values(row_index) = (*this)(xpos, ypos, zpos);
		//}

		std::cerr << "finished evaluation" << std::endl;
	}

	template<Degree deg>
	void TensorBSplines<deg>::precompute_pi2(Matrix &pi2) {
		const Matrix &d2B = pfs2.B;
		pi2 = INTofMULT(d2B,d2B);
	}

	template<Degree deg>
	void TensorBSplines<deg>::precompute_Pi2(Matrix &Pi2, const Size N) {
		Matrix pi2;
		precompute_pi2(pi2);

		//std::cerr << "pi2 = \n" << pi2 << std::endl;

		Pi2 = Matrix::Zero(N, N);
		Scalar delta = (Scalar) 1.0/(N-deg);
		for (Index i = 0; i < N-deg; i++) {
			for (Index r = 0; r <= deg; r++) {
				for (Index s = 0; s <= deg; s++) {
					Pi2(i+r, i+s) += pi2(r, s) * delta;
				}
			}
		}
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_smooth_1d_mat(SparseMatrix &smooth_1d_mat, const Size N) {

		std::cerr << "\nprepare smooth_1d_mat:" << std::endl;

		Matrix Pi2;
		precompute_Pi2(Pi2, N);

		//std::cerr << "Pi2 = \n" << Pi2 << std::endl;

		std::vector<Triplet> smooth_tripleList;
		for (Index row_index = 0; row_index < N; row_index++) {
			for (Index column_index = std::max(row_index-deg, 0); column_index <= std::min(row_index+deg, N-1); column_index++) {
				Scalar weight = Pi2(row_index, column_index);
				smooth_tripleList.push_back(Triplet(row_index, column_index, weight));
			}
		}
		smooth_1d_mat.resize(N, N);
		smooth_1d_mat.setFromTriplets(smooth_tripleList.begin(), smooth_tripleList.end());

		std::cerr << "finished smooth_mat" << std::endl;
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_bs_iws(IndexWeightVec &iws, Scalar pos, const Size N) {
		Scalar delta = (Scalar) 1.0/(N-deg);

		Index i = (Index)(pos/delta);
		i = (i>=N-deg) ? (N-deg-1) : i;

		const Matrix &B = pfs0.B;
		Scalar u = pos/delta - i;
		Vector pu = Vector::Ones(deg+1);
		for (Index h = 0; h < deg; h++) pu(h+1) = pu(h) * u;
		Vector bu;
		bu = B*pu;

		for (Index r = 0; r <= deg; r++) {
			Index index = i+r;
			iws.push_back(std::make_pair(index, bu(r)));
		}
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_data_1d_mat(const Vector &points_1d, SparseMatrix &data_1d_mat, const Size N) {

		std::cerr << "\nprepare data_1d_mat:" << std::endl;

		Size npts = (Size) points_1d.size();
		data_1d_mat.resize(npts, N);

		std::vector<Triplet> data_list;
		for (Index row_index = 0; row_index < npts; row_index++) {
			Scalar pos = points_1d(row_index);

			IndexWeightVec iws;
			make_bs_iws(iws, pos, N);
			for(auto it = iws.begin(); it!= iws.end(); it++) {
				Index column_index = it->first;
				Scalar weight = it->second;
				data_list.push_back(Triplet(row_index, column_index, weight));
			}
		}
		data_1d_mat.setFromTriplets(data_list.begin(), data_list.end());

		std::cerr << "finished data_1d_mat" << std::endl;
	}

	template<Degree deg>
	void TensorBSplines<deg>::precompute_pi0pi1pi2(Matrix &pi0, Matrix &pi1, Matrix &pi2) {
		const Matrix &B = pfs0.B, &dB = pfs1.B, &d2B = pfs2.B;
		pi0 = INTofMULT(B,B);
		pi1 = INTofMULT(dB,dB);
		pi2 = INTofMULT(d2B,d2B);
	}

	template<Degree deg>
	void TensorBSplines<deg>::precompute_Pi0Pi1Pi2(Matrix &Pi0, Matrix &Pi1, Matrix &Pi2, const Size N) {
		Matrix pi0, pi1, pi2;
		precompute_pi0pi1pi2(pi0, pi1, pi2);

		//std::cerr << "pi0 = \n" << pi0 << std::endl;
		//std::cerr << "pi1 = \n" << pi1 << std::endl;
		//std::cerr << "pi2 = \n" << pi2 << std::endl;

		Pi0 = Matrix::Zero(N, N), Pi1 = Matrix::Zero(N, N), Pi2 = Matrix::Zero(N, N);
		Scalar delta = (Scalar) 1.0/(N-deg);
		for (Index i = 0; i < N-deg; i++) {
			for (Index r = 0; r <= deg; r++) {
				for (Index s = 0; s <= deg; s++) {
					Pi0(i+r, i+s) += pi0(r, s) * delta;
					Pi1(i+r, i+s) += pi1(r, s) * delta;
					Pi2(i+r, i+s) += pi2(r, s) * delta;
				}
			}
		}
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_smooth_mat(SparseMatrix &smooth_mat, const Size N) {

		std::cerr << "\nprepare smooth_mat:" << std::endl;

		Matrix Pi0, Pi1, Pi2;
		precompute_Pi0Pi1Pi2(Pi0, Pi1, Pi2, N);

		//std::cerr << "Pi0 = \n" << Pi0 << std::endl;
		//std::cerr << "Pi1 = \n" << Pi1 << std::endl;
		//std::cerr << "Pi2 = \n" << Pi2 << std::endl;

		std::vector<Triplet> smooth_tripleList;
		for (Index k = 0; k < N; k++) {
			for (Index j = 0; j < N; j++) {
				for (Index i = 0; i < N; i++) {
					Index row_index = k*N*N + j*N + i;
					for (Index kk = std::max(k-deg, 0); kk <= std::min(k+deg, N-1); kk++) {
						for (Index jj = std::max(j-deg, 0); jj <= std::min(j+deg, N-1); jj++) {
							for (Index ii = std::max(i-deg, 0); ii <= std::min(i+deg, N-1); ii++) {
								Index column_index = kk*N*N + jj*N + ii;
								if (column_index < row_index) continue;
								Scalar weight = Pi2(i, ii) * Pi0(j, jj) * Pi0(k, kk) 
									+ Pi0(i, ii) * Pi2(j, jj) * Pi0(k, kk)
									+ Pi0(i, ii) * Pi0(j, jj) * Pi2(k, kk)
									+ Pi1(i, ii) * Pi1(j, jj) * Pi0(k, kk) * 2
									+ Pi1(i, ii) * Pi0(j, jj) * Pi1(k, kk) * 2
									+ Pi0(i, ii) * Pi1(j, jj) * Pi1(k, kk) * 2;
								smooth_tripleList.push_back(Triplet(row_index, column_index, weight));
								if (column_index > row_index) smooth_tripleList.push_back(Triplet(column_index, row_index, weight));
							}
						}
					}
				}
			}
		}
		smooth_mat.resize(N*N*N, N*N*N);
		smooth_mat.setFromTriplets(smooth_tripleList.begin(), smooth_tripleList.end());

		std::cerr << "finished smooth_mat" << std::endl;
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_tbs_iws(IndexWeightVec &iws, Scalar xpos, Scalar ypos, Scalar zpos, const Size N) {
		Scalar delta = (Scalar) 1.0/(N-deg);

		Index i = (Index)(xpos/delta), j = (Index)(ypos/delta), k = (Index)(zpos/delta);
		i = (i>=N-deg) ? (N-deg-1) : i;
		j = (j>=N-deg) ? (N-deg-1) : j;
		k = (k>=N-deg) ? (N-deg-1) : k;
		Scalar u = xpos/delta - i, v = ypos/delta - j, w = zpos/delta - k;

		const Matrix &B = pfs0.B;
		Vector pu = Vector::Ones(deg+1), pv = Vector::Ones(deg+1), pw = Vector::Ones(deg+1);
		for(Index h = 0; h < deg; h++) {
			pu(h+1) = pu(h) * u;
			pv(h+1) = pv(h) * v;
			pw(h+1) = pw(h) * w;
		}
		Vector bu, bv, bw;
		bu = B*pu;
		bv = B*pv;
		bw = B*pw;

		for(Index t = 0; t <= deg; t++) {
			for(Index s = 0; s <= deg; s++) {
				for(Index r = 0; r <= deg; r++) {
					Index index = i+r+(j+s)*N+(k+t)*N*N;
					iws.push_back(std::make_pair(index, bu(r)*bv(s)*bw(t)));
				}
			}
		}
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_data_mat(const PointSet &points, SparseMatrix &data_mat, const Size N) {

		std::cerr << "\nprepare data_mat:" << std::endl;

		Size npts = (Size) points.rows();
		data_mat.resize(npts, N*N*N);

		std::vector<Triplet> data_list;
		Scalar delta = (Scalar) 1.0/(N-deg);
		Arrayi ijk = (points/delta).array().floor().cast<Index>();
		for(Index h = 0; h < npts; h++) ijk(h,0) = ijk(h,0) < (N-deg) ? ijk(h,0) : (N-deg-1);
		for(Index h = 0; h < npts; h++) ijk(h,1) = ijk(h,1) < (N-deg) ? ijk(h,1) : (N-deg-1);
		for(Index h = 0; h < npts; h++) ijk(h,2) = ijk(h,2) < (N-deg) ? ijk(h,2) : (N-deg-1);

		Array uvw = (points/delta).array() - ijk.cast<Scalar>();

		const Matrix &B = pfs0.B;
		Array U = Array::Ones(npts, deg+1), V = Array::Ones(npts, deg+1), W = Array::Ones(npts, deg+1);
		for (Index h = 0; h < deg; h++) {
			U.col(h+1) = U.col(h) * uvw.col(0);
			V.col(h+1) = V.col(h) * uvw.col(1);
			W.col(h+1) = W.col(h) * uvw.col(2);
		}

		Matrix bu = U.matrix()*B.transpose(), bv = V.matrix()*B.transpose(), bw = W.matrix()*B.transpose();

		for (Index kk = 0; kk <= deg; kk++) {
			for (Index jj = 0; jj <= deg; jj++) {
				for (Index ii = 0; ii <= deg; ii++) {
					Arrayi indices = N*N*(ijk.col(2)+kk) + N*(ijk.col(1)+jj) + ijk.col(0)+ii;
					Array bubvbw = bu.col(ii).array()*bv.col(jj).array()*bw.col(kk).array();
					for(Index h = 0; h < npts; h++) {
						Index row_index = h, column_index = indices(h);
						Scalar weight = bubvbw(h);
						data_list.push_back(Triplet(row_index, column_index, weight));
					}
				}
			}
		}
		//for (Index row_index = 0; row_index < npts; row_index++) {
		//	Point point = points.row(row_index);
		//	Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();

		//	IndexWeightVec iws;
		//	make_tbs_iws(iws, xpos, ypos, zpos);
		//	for(auto it = iws.begin(); it!= iws.end(); it++) {
		//		Index column_index = it->first;
		//		Scalar weight = it->second;
		//		data_list.push_back(Triplet(row_index, column_index, weight));
		//	}
		//}
		data_mat.setFromTriplets(data_list.begin(), data_list.end());

		std::cerr << "finished data_mat" << std::endl;
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_tbs_iws_duvw_ws(IndexWeightVec &iws, std::vector<Scalar> &du_ws,
		std::vector<Scalar> &dv_ws, std::vector<Scalar> &dw_ws, Scalar xpos, Scalar ypos, Scalar zpos, const Size N) {
			Scalar delta = (Scalar) 1.0/(N-deg);

			Index i = (Index)(xpos/delta), j = (Index)(ypos/delta), k = (Index)(zpos/delta);
			i = (i>=N-deg) ? (N-deg-1) : i;
			j = (j>=N-deg) ? (N-deg-1) : j;
			k = (k>=N-deg) ? (N-deg-1) : k;
			Scalar u = xpos/delta - i, v = ypos/delta - j, w = zpos/delta - k;

			const Matrix &B = pfs0.B;
			const Matrix &dB = pfs1.B;
			Vector pu = Vector::Ones(deg+1), pv = Vector::Ones(deg+1), pw = Vector::Ones(deg+1);
			for(Index h = 0; h < deg; h++) {
				pu(h+1) = pu(h) * u;
				pv(h+1) = pv(h) * v;
				pw(h+1) = pw(h) * w;
			}

			Vector bu, bv, bw, dbu, dbv, dbw;
			bu = B*pu;
			bv = B*pv;
			bw = B*pw;
			dbu = dB*pu;
			dbv = dB*pv;
			dbw = dB*pw;

			for(Index t = 0; t <= deg; t++) {
				for(Index s = 0; s <= deg; s++) {
					for(Index r = 0; r <= deg; r++) {
						Index index = i+r+(j+s)*N+(k+t)*N*N;
						iws.push_back(std::make_pair(index, bu(r)*bv(s)*bw(t)));
						du_ws.push_back(dbu(r)*bv(s)*bw(t));
						dv_ws.push_back(bu(r)*dbv(s)*bw(t));
						dw_ws.push_back(bu(r)*bv(s)*dbw(t));
					}
				}
			}
	}

	template<Degree deg>
	void TensorBSplines<deg>::make_data_duvw_mat(const PointSet &points, SparseMatrix &data_mat, SparseMatrix &du_mat, 
		SparseMatrix &dv_mat, SparseMatrix &dw_mat, const Size N) {

			std::cerr << "prepare data_mat, du_mat, dv_mat, dw_mat" << std::endl;

			Size npts = (Size) points.rows();
			data_mat.resize(npts, N*N*N);
			du_mat.resize(npts, N*N*N);
			dv_mat.resize(npts, N*N*N);
			dw_mat.resize(npts, N*N*N);

			std::vector<Triplet> data_list, du_list, dv_list, dw_list;
			Scalar delta = (Scalar) 1.0/(N-deg);
			Arrayi ijk = (points/delta).array().floor().cast<Index>();
			for(Index h = 0; h < npts; h++) ijk(h,0) = ijk(h,0) < (N-deg) ? ijk(h,0) : (N-deg-1);
			for(Index h = 0; h < npts; h++) ijk(h,1) = ijk(h,1) < (N-deg) ? ijk(h,1) : (N-deg-1);
			for(Index h = 0; h < npts; h++) ijk(h,2) = ijk(h,2) < (N-deg) ? ijk(h,2) : (N-deg-1);

			Array uvw = (points/delta).array() - ijk.cast<Scalar>();

			const Matrix &B = pfs0.B;
			const Matrix &dB = pfs1.B;
			Array U = Array::Ones(npts, deg+1), V = Array::Ones(npts, deg+1), W = Array::Ones(npts, deg+1);
			for (Index h = 0; h < deg; h++) {
				U.col(h+1) = U.col(h) * uvw.col(0);
				V.col(h+1) = V.col(h) * uvw.col(1);
				W.col(h+1) = W.col(h) * uvw.col(2);
			}

			Matrix bu = U.matrix()*B.transpose(), bv = V.matrix()*B.transpose(), bw = W.matrix()*B.transpose();
			Matrix dbu = U.matrix()*dB.transpose(), dbv = V.matrix()*dB.transpose(), dbw = W.matrix()*dB.transpose();

			for (Index kk = 0; kk <= deg; kk++) {
				for (Index jj = 0; jj <= deg; jj++) {
					for (Index ii = 0; ii <= deg; ii++) {
						Arrayi indices = N*N*(ijk.col(2)+kk) + N*(ijk.col(1)+jj) + ijk.col(0)+ii;
						Array bubvbw = bu.col(ii).array()*bv.col(jj).array()*bw.col(kk).array();
						Array dbubvbw = dbu.col(ii).array()*bv.col(jj).array()*bw.col(kk).array();
						Array budbvbw = bu.col(ii).array()*dbv.col(jj).array()*bw.col(kk).array();
						Array bubvdbw = bu.col(ii).array()*bv.col(jj).array()*dbw.col(kk).array();
						for(Index h = 0; h < npts; h++) {
							Index row_index = h, column_index = indices(h);
							Scalar weight = bubvbw(h), du_weight = dbubvbw(h), dv_weight = budbvbw(h), dw_weight = bubvdbw(h);
							data_list.push_back(Triplet(row_index, column_index, weight));
							du_list.push_back(Triplet(row_index, column_index, du_weight));
							dv_list.push_back(Triplet(row_index, column_index, dv_weight));
							dw_list.push_back(Triplet(row_index, column_index, dw_weight));
						}
					}
				}
			}
			//for (Index row_index = 0; row_index < npts; row_index++) {
			//	Point point = points.row(row_index);
			//	Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();

			//	IndexWeightVec iws;
			//	std::vector<Scalar> du_ws, dv_ws, dw_ws;
			//	make_tbs_iws_duvw_ws(iws, du_ws, dv_ws, dw_ws, xpos, ypos, zpos, N);
			//	Index h = 0;
			//	for(auto it = iws.begin(); it!= iws.end(); it++, h++) {
			//		Index column_index = it->first;
			//		Scalar weight = it->second, du_weight = du_ws[h], dv_weight = dv_ws[h], dw_weight = dw_ws[h];
			//		data_list.push_back(Triplet(row_index, column_index, weight));
			//		du_list.push_back(Triplet(row_index, column_index, du_weight));
			//		dv_list.push_back(Triplet(row_index, column_index, dv_weight));
			//		dw_list.push_back(Triplet(row_index, column_index, dw_weight));
			//	}
			//}
			data_mat.setFromTriplets(data_list.begin(), data_list.end());
			du_mat.setFromTriplets(du_list.begin(), du_list.end());
			dv_mat.setFromTriplets(dv_list.begin(), dv_list.end());
			dw_mat.setFromTriplets(dw_list.begin(), dw_list.end());

			std::cerr << "finished data_mat, du_mat, dv_mat, dw_mat" << std::endl;
	}
};
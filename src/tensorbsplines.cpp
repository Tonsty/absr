#include <tensorbsplines.h>
#include <iostream>

using namespace absr;

Size N = 128;

static struct NormalizedCubicSpline {
	Matrix B, dB, d2B;
	NormalizedCubicSpline() {
		B.resize(4, 4);
		dB.resize(4, 4);
		d2B.resize(4, 4);
		B<< 1, -3, 3, -1,
			4, 0, -6, 3,
			1, 3, 3, -3,
			0, 0, 0, 1;
		dB<< -3, 6, -3, 0,
			0, -12, 9, 0,
			3, 6, -9, 0,
			0, 0, 3, 0;
		d2B<< 6, -6, 0, 0,
			-12, 18, 0, 0,
			6, -18, 0, 0,
			0, 6, 0, 0;
		B/=6; dB/=6; d2B/=6;	
	}
} ncs;

Vector MULT(const Vector &U, const Vector &V) {
	Vector W = Vector::Zero(7, 1);
	for(Index i = 0; i < 4; i++) {
		for (Index j = 0; j < 4; j++) {
			W(i+j) += U(i)*V(j);
		}
	}
	return W;
}

Matrix INTofMULT(const Matrix &A, const Matrix &B) {
	Matrix C = Matrix::Zero(4, 4);
	for (Index i = 0; i < 4; i++) {
		for (Index j = 0; j < 4; j++) {
			Vector W = MULT(A.row(i), B.row(j));
			//std::cerr << "W=\n" << W << std::endl; 
			Scalar integral_sum = 0;
			for(Index k = 0; k < 7; k++) {
				integral_sum += W(k)/(k+1);
			}
			C(i,j) = integral_sum;
		}
	}
	return C;
}

void precompute_pi2(Matrix &pi2) {
	const Matrix &d2B = ncs.d2B;
	pi2 = INTofMULT(d2B,d2B);
}

void precompute_Pi2(Matrix &Pi2) {
	Matrix pi2;
	precompute_pi2(pi2);

	//std::cerr << "pi2 = \n" << pi2 << std::endl;

	Pi2 = Matrix::Zero(N, N);
	Scalar delta = (Scalar) 1.0/(N-3);
	for (Index i = 0; i < N-3; i++) {
		for (Index r = 0; r <= 3; r++) {
			for (Index s = 0; s <= 3; s++) {
				Pi2(i+r, i+s) += pi2(r, s) * delta;
			}
		}
	}
}

void TensorBSplines::make_smooth_1d_mat(SparseMatrix &smooth_1d_mat) {

	std::cerr << "\nprepare smooth_1d_mat:" << std::endl;

	Matrix Pi2;
	precompute_Pi2(Pi2);

	//std::cerr << "Pi2 = \n" << Pi2 << std::endl;

	std::vector<Triplet> smooth_tripleList;
	for (Index row_index = 0; row_index < N; row_index++) {
		for (Index column_index = std::max(row_index-3, 0); column_index <= std::min(row_index+3, N-1); column_index++) {
			Scalar weight = Pi2(row_index, column_index);
			smooth_tripleList.push_back(Triplet(row_index, column_index, weight));
		}
	}
	smooth_1d_mat.resize(N, N);
	smooth_1d_mat.setFromTriplets(smooth_tripleList.begin(), smooth_tripleList.end());

	std::cerr << "finished smooth_mat" << std::endl;
}

void make_bs_iws(std::vector<std::pair<Index, Scalar>> &iws, Scalar pos) {
	Scalar delta = (Scalar) 1.0/(N-3);

	Index i = (Index)(pos/delta);
	i = (i>=N-3) ? (N-4) : i;

	const Matrix &B = ncs.B;
	Scalar u = pos/delta - i;
	Vector pu(4, 1);
	pu << 1, u, u*u, u*u*u;
	Vector bu;
	bu = B*pu;

	for (Index r = 0; r <= 3; r++) {
		Index index = i+r;
		iws.push_back(std::make_pair(index, bu(r)));
	}
}

void TensorBSplines::make_data_1d_mat(const Vector &points_1d, SparseMatrix &data_1d_mat) {

	std::cerr << "\nprepare data_1d_mat:" << std::endl;

	Size npts = (Size) points_1d.size();
	data_1d_mat.resize(npts, N);

	std::vector<Triplet> data_list;
	for (Index row_index = 0; row_index < npts; row_index++) {
		Scalar pos = points_1d(row_index);

		std::vector<std::pair<Index, Scalar>> iws;
		make_bs_iws(iws, pos);
		for(auto it = iws.begin(); it!= iws.end(); it++) {
			Index column_index = it->first;
			Scalar weight = it->second;
			data_list.push_back(Triplet(row_index, column_index, weight));
		}
	}
	data_1d_mat.setFromTriplets(data_list.begin(), data_list.end());

	std::cerr << "finished data_1d_mat" << std::endl;
}

void precompute_pi0pi1pi2(Matrix &pi0, Matrix &pi1, Matrix &pi2) {
	const Matrix &B = ncs.B, &dB = ncs.dB, &d2B = ncs.d2B;
	pi0 = INTofMULT(B,B);
	pi1 = INTofMULT(dB,dB);
	pi2 = INTofMULT(d2B,d2B);
}

void precompute_Pi0Pi1Pi2(Matrix &Pi0, Matrix &Pi1, Matrix &Pi2) {
	Matrix pi0, pi1, pi2;
	precompute_pi0pi1pi2(pi0, pi1, pi2);

	//std::cerr << "pi0 = \n" << pi0 << std::endl;
	//std::cerr << "pi1 = \n" << pi1 << std::endl;
	//std::cerr << "pi2 = \n" << pi2 << std::endl;

	Pi0 = Matrix::Zero(N, N), Pi1 = Matrix::Zero(N, N), Pi2 = Matrix::Zero(N, N);
	Scalar delta = (Scalar) 1.0/(N-3);
	for (Index i = 0; i < N-3; i++) {
		for (Index r = 0; r <= 3; r++) {
			for (Index s = 0; s <= 3; s++) {
				Pi0(i+r, i+s) += pi0(r, s) * delta;
				Pi1(i+r, i+s) += pi1(r, s) * delta;
				Pi2(i+r, i+s) += pi2(r, s) * delta;
			}
		}
	}
}

void TensorBSplines::make_smooth_mat(SparseMatrix &smooth_mat) {

	std::cerr << "\nprepare smooth_mat:" << std::endl;

	Matrix Pi0, Pi1, Pi2;
	precompute_Pi0Pi1Pi2(Pi0, Pi1, Pi2);

	//std::cerr << "Pi0 = \n" << Pi0 << std::endl;
	//std::cerr << "Pi1 = \n" << Pi1 << std::endl;
	//std::cerr << "Pi2 = \n" << Pi2 << std::endl;

	std::vector<Triplet> smooth_tripleList;
	for (Index k = 0; k < N; k++) {
		for (Index j = 0; j < N; j++) {
			for (Index i = 0; i < N; i++) {
				Index row_index = k*N*N + j*N + i;
				for (Index kk = std::max(k-3, 0); kk <= std::min(k+3, N-1); kk++) {
					for (Index jj = std::max(j-3, 0); jj <= std::min(j+3, N-1); jj++) {
						for (Index ii = std::max(i-3, 0); ii <= std::min(i+3, N-1); ii++) {
							Index column_index = kk*N*N + jj*N + ii;
							Scalar weight = Pi2(i, ii) * Pi0(j, jj) * Pi0(k, kk) 
								+ Pi0(i, ii) * Pi2(j, jj) * Pi0(k, kk)
								+ Pi0(i, ii) * Pi0(j, jj) * Pi2(k, kk)
								+ Pi1(i, ii) * Pi1(j, jj) * Pi0(k, kk) * 2
								+ Pi1(i, ii) * Pi0(j, jj) * Pi1(k, kk) * 2
								+ Pi0(i, ii) * Pi1(j, jj) * Pi1(k, kk) * 2;
							smooth_tripleList.push_back(Triplet(row_index, column_index, weight));
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

void make_tbs_iws(std::vector<std::pair<Index, Scalar>> &iws, Scalar xpos, Scalar ypos, Scalar zpos) {
	Scalar delta = (Scalar) 1.0/(N-3);

	Index i = (Index)(xpos/delta), j = (Index)(ypos/delta), k = (Index)(zpos/delta);
	i = (i>=N-3) ? (N-4) : i;
	j = (j>=N-3) ? (N-4) : j;
	k = (k>=N-3) ? (N-4) : k;
	Scalar u = xpos/delta - i, v = ypos/delta - j, w = zpos/delta - k;

	const Matrix &B = ncs.B;
	Vector pu(4, 1), pv(4, 1), pw(4, 1);
	pu << 1, u, u*u, u*u*u;
	pv << 1, v, v*v, v*v*v;
	pw << 1, w, w*w, w*w*w;
	Vector bu, bv, bw;
	bu = B*pu;
	bv = B*pv;
	bw = B*pw;

	for(Index t = 0; t <= 3; t++) {
		for(Index s = 0; s <= 3; s++) {
			for(Index r = 0; r <= 3; r++) {
				Index index = i+r+(j+s)*N+(k+t)*N*N;
				iws.push_back(std::make_pair(index, bu(r)*bv(s)*bw(t)));
			}
		}
	}
}

void TensorBSplines::make_data_mat(const PointSet &points, SparseMatrix &data_mat) {

	std::cerr << "\nprepare data_mat:" << std::endl;

	Size npts = (Size) points.rows();
	data_mat.resize(npts, N*N*N);

	std::vector<Triplet> data_list;
	Scalar delta = (Scalar) 1.0/(N-3);
	Arrayi ijk = (points/delta).array().floor().cast<Index>();
	for(Index h = 0; h < npts; h++) ijk(h,0) = ijk(h,0) < (N-3) ? ijk(h,0) : (N-4);
	for(Index h = 0; h < npts; h++) ijk(h,1) = ijk(h,1) < (N-3) ? ijk(h,1) : (N-4);
	for(Index h = 0; h < npts; h++) ijk(h,2) = ijk(h,2) < (N-3) ? ijk(h,2) : (N-4);

	Array uvw = (points/delta).array() - ijk.cast<Scalar>();

	const Matrix &B = ncs.B;
	Array U(npts, 4), V(npts, 4), W(npts, 4);
	U.col(0) = V.col(0) = W.col(0) = Array::Ones(npts, 1);
	U.col(1) = uvw.col(0);        V.col(1) = uvw.col(1);        W.col(1) = uvw.col(2);
	U.col(2) = U.col(1)*U.col(1); V.col(2) = V.col(1)*V.col(1); W.col(2) = W.col(1)*W.col(1);
	U.col(3) = U.col(1)*U.col(2); V.col(3) = V.col(1)*V.col(2); W.col(3) = W.col(1)*W.col(2);

	Matrix bu = U.matrix()*B.transpose(), bv = V.matrix()*B.transpose(), bw = W.matrix()*B.transpose();

	for (Index kk = 0; kk <= 3; kk++) {
		for (Index jj = 0; jj <= 3; jj++) {
			for (Index ii = 0; ii <= 3; ii++) {
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

	//	std::vector<std::pair<Index, Scalar>> iws;
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

void make_tbs_iws_duvw_ws(std::vector<std::pair<Index, Scalar>> &iws, 
						  std::vector<Scalar> &du_ws,
						  std::vector<Scalar> &dv_ws,
						  std::vector<Scalar> &dw_ws,
						  Scalar xpos, Scalar ypos, Scalar zpos) {
	Scalar delta = (Scalar) 1.0/(N-3);

	Index i = (Index)(xpos/delta), j = (Index)(ypos/delta), k = (Index)(zpos/delta);
	i = (i>=N-3) ? (N-4) : i;
	j = (j>=N-3) ? (N-4) : j;
	k = (k>=N-3) ? (N-4) : k;
	Scalar u = xpos/delta - i, v = ypos/delta - j, w = zpos/delta - k;

	const Matrix &B = ncs.B;
	const Matrix &dB = ncs.dB;
	Vector pu(4, 1), pv(4, 1), pw(4, 1);
	pu << 1, u, u*u, u*u*u;
	pv << 1, v, v*v, v*v*v;
	pw << 1, w, w*w, w*w*w;

	Vector bu, bv, bw, dbu, dbv, dbw;
	bu = B*pu;
	bv = B*pv;
	bw = B*pw;
	dbu = dB*pu;
	dbv = dB*pv;
	dbw = dB*pw;

	for(Index t = 0; t <= 3; t++) {
		for(Index s = 0; s <= 3; s++) {
			for(Index r = 0; r <= 3; r++) {
				Index index = i+r+(j+s)*N+(k+t)*N*N;
				iws.push_back(std::make_pair(index, bu(r)*bv(s)*bw(t)));
				du_ws.push_back(dbu(r)*bv(s)*bw(t));
				dv_ws.push_back(bu(r)*dbv(s)*bw(t));
				dw_ws.push_back(bu(r)*bv(s)*dbw(t));
			}
		}
	}
}

void TensorBSplines::make_data_duvw_mat(const PointSet &points, 
						SparseMatrix &data_mat, 
						SparseMatrix &du_mat, 
						SparseMatrix &dv_mat, 
						SparseMatrix &dw_mat) {

	std::cerr << "prepare data_mat, du_mat, dv_mat, dw_mat" << std::endl;

	Size npts = (Size) points.rows();
	data_mat.resize(npts, N*N*N);
	du_mat.resize(npts, N*N*N);
	dv_mat.resize(npts, N*N*N);
	dw_mat.resize(npts, N*N*N);

	std::vector<Triplet> data_list, du_list, dv_list, dw_list;
	Scalar delta = (Scalar) 1.0/(N-3);
	Arrayi ijk = (points/delta).array().floor().cast<Index>();
	for(Index h = 0; h < npts; h++) ijk(h,0) = ijk(h,0) < (N-3) ? ijk(h,0) : (N-4);
	for(Index h = 0; h < npts; h++) ijk(h,1) = ijk(h,1) < (N-3) ? ijk(h,1) : (N-4);
	for(Index h = 0; h < npts; h++) ijk(h,2) = ijk(h,2) < (N-3) ? ijk(h,2) : (N-4);

	Array uvw = (points/delta).array() - ijk.cast<Scalar>();

	const Matrix &B = ncs.B;
	const Matrix &dB = ncs.dB;
	Array U(npts, 4), V(npts, 4), W(npts, 4);
	U.col(0) = V.col(0) = W.col(0) = Array::Ones(npts, 1);
	U.col(1) = uvw.col(0);        V.col(1) = uvw.col(1);        W.col(1) = uvw.col(2);
	U.col(2) = U.col(1)*U.col(1); V.col(2) = V.col(1)*V.col(1); W.col(2) = W.col(1)*W.col(1);
	U.col(3) = U.col(1)*U.col(2); V.col(3) = V.col(1)*V.col(2); W.col(3) = W.col(1)*W.col(2);

	Matrix bu = U.matrix()*B.transpose(), bv = V.matrix()*B.transpose(), bw = W.matrix()*B.transpose();
	Matrix dbu = U.matrix()*dB.transpose(), dbv = V.matrix()*dB.transpose(), dbw = W.matrix()*dB.transpose();

	for (Index kk = 0; kk <= 3; kk++) {
		for (Index jj = 0; jj <= 3; jj++) {
			for (Index ii = 0; ii <= 3; ii++) {
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

	//	std::vector<std::pair<Index, Scalar>> iws;
	//	std::vector<Scalar> du_ws, dv_ws, dw_ws;
	//	make_tbs_iws_duvw_ws(iws, du_ws, dv_ws, dw_ws, xpos, ypos, zpos);
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


void TensorBSplines::evaluate(const PointSet &points, Vector &values) {
	Size npts = (Size) points.rows();
	values = Vector::Zero(npts);

	std::cerr << "begin evaluation:" << std::endl;

	Scalar delta = (Scalar) 1.0/(N-3);
	Arrayi ijk = (points/delta).array().floor().cast<Index>();
	for(Index h = 0; h < npts; h++) ijk(h,0) = ijk(h,0) < (N-3) ? ijk(h,0) : (N-4);
	for(Index h = 0; h < npts; h++) ijk(h,1) = ijk(h,1) < (N-3) ? ijk(h,1) : (N-4);
	for(Index h = 0; h < npts; h++) ijk(h,2) = ijk(h,2) < (N-3) ? ijk(h,2) : (N-4);

	Array uvw = (points/delta).array() - ijk.cast<Scalar>();

	const Matrix &B = ncs.B;
	Array U(npts, 4), V(npts, 4), W(npts, 4);
	U.col(0) = V.col(0) = W.col(0) = Array::Ones(npts, 1);
	U.col(1) = uvw.col(0);        V.col(1) = uvw.col(1);        W.col(1) = uvw.col(2);
	U.col(2) = U.col(1)*U.col(1); V.col(2) = V.col(1)*V.col(1); W.col(2) = W.col(1)*W.col(1);
	U.col(3) = U.col(1)*U.col(2); V.col(3) = V.col(1)*V.col(2); W.col(3) = W.col(1)*W.col(2);

	Matrix bu = U.matrix()*B.transpose(), bv = V.matrix()*B.transpose(), bw = W.matrix()*B.transpose();

	for (Index kk = 0; kk <= 3; kk++) {
		for (Index jj = 0; jj <= 3; jj++) {
			for (Index ii = 0; ii <= 3; ii++) {
				Arrayi indices = N*N*(ijk.col(2)+kk) + N*(ijk.col(1)+jj) + ijk.col(0)+ii;
				Array bubvbw = bu.col(ii).array()*bv.col(jj).array()*bw.col(kk).array();
				for(Index h = 0; h < npts; h++) values(h) += controls_(indices(h))*bubvbw(h);
			}
		}
	}

	//for (Index row_index = 0; row_index < npts; row_index++) {
	//	Point point = points.row(row_index);
	//	Scalar xpos = point.x(), ypos = point.y(), zpos = point.z();

	//	std::vector<std::pair<Index, Scalar>> iws;
	//	make_tbs_iws(iws, xpos, ypos, zpos);

	//	Scalar sum_value = 0;
	//	for(auto it = iws.begin(); it!= iws.end(); it++) {
	//		Index control_index = it->first;
	//		Scalar weight = it->second;
	//		sum_value += controls_(control_index) * weight;
	//	}
	//	values(row_index) = sum_value;
	//}

	std::cerr << "finished evaluation" << std::endl;
}
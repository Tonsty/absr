#include <activetbsfitting.h>
#include <iostream>

using namespace absr;

void ActiveTBSFitting3L::fitting(ActiveTBS &atbs) {
	PointSet &points = points_;
	NormalSet &normals = normals_;
	const Scalar lambda = lambda_;
	const Scalar epsilon = epsilon_;
	const Size N = atbs.N_;

	Size npts = (Size) points.rows();
	PointSet points_3L(npts*3, 3);
	points_3L << points, points+epsilon*normals, points-epsilon*normals;
	Vector values_3L(npts*3);
	values_3L << Vector::Zero(npts), Vector::Ones(npts), -Vector::Ones(npts);

	ActiveTBS::generate_active_map(atbs.amp_, points_3L, N);

	SparseMatrix data_mat;
	ActiveTBS::make_data_mat(atbs.amp_, points_3L, data_mat, N);

	SparseMatrix smooth_mat; 
	ActiveTBS::make_smooth_mat(atbs.amp_, smooth_mat, N);

	SparseMatrix left_mat = data_mat.transpose() * data_mat + smooth_mat * lambda * npts * 3;
	Vector right_vec = data_mat.transpose() * values_3L;

	std::cerr << "\nbegin solving active 3L:" << std::endl;
	std::cerr << "equation number = 1, equation size = " << 3*npts << " * " << atbs.amp_.size() << std::endl;

	Eigen::ConjugateGradient<SparseMatrix> cg;
	cg.compute(left_mat);
	atbs.controls_ = cg.solve(right_vec);

	std::cerr << "finished solving active 3L" << std::endl;
}

void ActiveTBSFittingJuttler::fitting(ActiveTBS &atbs) {
	PointSet &points = points_;
	NormalSet &normals = normals_;
	const Scalar lambda = lambda_;
	const Scalar kappa = kappa_;
	const Size N = atbs.N_;

	ActiveTBS::generate_active_map(atbs.amp_, points, N);

	SparseMatrix data_mat, du_mat, dv_mat, dw_mat;
	ActiveTBS::make_data_duvw_mat(atbs.amp_, points, data_mat, du_mat, dv_mat, dw_mat, N);

	Size npts = (Size) points.rows();
	SparseMatrix smooth_mat; 
	ActiveTBS::make_smooth_mat(atbs.amp_, smooth_mat, N);

	SparseMatrix left_mat = data_mat.transpose() * data_mat + 
		(du_mat.transpose() * du_mat + 
		dv_mat.transpose() * dv_mat + 
		dw_mat.transpose() * dw_mat) * kappa +
		smooth_mat * lambda * npts;
	Vector right_vec = (du_mat.transpose() * normals.col(0) +
		dv_mat.transpose() * normals.col(1) +
		dw_mat.transpose() * normals.col(2)) * kappa;

	std::cerr << "\nbegin solving active Juttler" << std::endl;
	std::cerr << "equation number = 1, equation size = " << npts << " * " << atbs.amp_.size() << std::endl;

	Eigen::ConjugateGradient<SparseMatrix> cg;
	cg.compute(left_mat);
	atbs.controls_ = cg.solve(right_vec);

	std::cerr << "finished solving active Juttler" << std::endl;
}


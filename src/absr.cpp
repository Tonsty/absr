#include <absr.h>
#include <iostream>

using namespace absr;

extern Size N;

void ABSRFitting3L::fitting(TensorBSplines &tbs) {
	PointSet &points = points_;
	NormalSet &normals = normals_;
	const Scalar &lambda = lambda_;
	const Scalar &epsilon = epsilon_;

	Size npts = (Size) points.rows();
	PointSet points_3L(npts*3, 3);
	points_3L << points, points+epsilon*normals, points-epsilon*normals;
	Vector values_3L(npts*3);
	values_3L << Vector::Zero(npts), Vector::Ones(npts), -Vector::Ones(npts);

	SparseMatrix data_mat;
	TensorBSplines::make_data_mat(points_3L, data_mat);

	SparseMatrix smooth_mat; 
	TensorBSplines::make_smooth_mat(smooth_mat);

	SparseMatrix left_mat = data_mat.transpose() * data_mat + smooth_mat * lambda * npts * 3;
	Vector right_vec = data_mat.transpose() * values_3L;

	std::cerr << "\nbegin solving 3L:" << std::endl;
	std::cerr << "equation number = 1, equation size = " << 3*npts << " * " << N*N*N << std::endl;

	Eigen::ConjugateGradient<SparseMatrix> cg;
	cg.compute(left_mat);
	tbs.controls_ = cg.solve(right_vec);

	std::cerr << "finished solving 3L" << std::endl;
}

void ABSRFittingJuttler::fitting(TensorBSplines &tbs) {
	PointSet &points = points_;
	NormalSet &normals = normals_;
	const Scalar &lambda = lambda_;
	const Scalar &kappa = kappa_;

	SparseMatrix data_mat, du_mat, dv_mat, dw_mat;
	TensorBSplines::make_data_duvw_mat(points, data_mat, du_mat, dv_mat, dw_mat);

	Size npts = (Size) points.rows();
	SparseMatrix smooth_mat; 
	TensorBSplines::make_smooth_mat(smooth_mat);

	SparseMatrix left_mat = data_mat.transpose() * data_mat + 
		(du_mat.transpose() * du_mat + 
		dv_mat.transpose() * dv_mat + 
		dw_mat.transpose() * dw_mat) * kappa +
		smooth_mat * lambda * npts;
	Vector right_vec = (du_mat.transpose() * normals.col(0) +
		dv_mat.transpose() * normals.col(1) +
		dw_mat.transpose() * normals.col(2)) * kappa;

	std::cerr << "\nbegin solving Juttler" << std::endl;
	std::cerr << "equation number = 1, equation size = " << npts << " * " << N*N*N << std::endl;

	Eigen::ConjugateGradient<SparseMatrix> cg;
	cg.compute(left_mat);
	tbs.controls_ = cg.solve(right_vec);

	std::cerr << "finished solving Juttler" << std::endl;
}

void ABSRFittingSDF::fitting(TensorBSplines &tbs) {
	const SDF &sdf = sdf_;
	const Scalar &lambda = lambda_;
	const bool &global_fitting = global_fitting_;

	const Vector &values = sdf.values_;
	const Size grid_size = sdf.grid_size_;
	const Scalar voxel_length = sdf.voxel_length_;

	if (global_fitting) {
		PointSet points;
		sdf.topoints(points);

		SparseMatrix data_mat;
		TensorBSplines::make_data_mat(points, data_mat);

		SparseMatrix smooth_mat; 
		TensorBSplines::make_smooth_mat(smooth_mat);

		SparseMatrix left_mat = data_mat.transpose() * data_mat + smooth_mat * lambda * grid_size * grid_size * grid_size;
		Vector right_vec = data_mat.transpose() * values;

		std::cerr << "\nbegin global solving:" << std::endl;
		std::cerr << "equation number = 1, equation size = " << grid_size*grid_size*grid_size << " * " << N*N*N << std::endl;

		Eigen::ConjugateGradient<SparseMatrix> cg;
		cg.compute(left_mat);
		tbs.controls_ = cg.solve(right_vec);

		std::cerr << "finished global solving" << std::endl;
	} else {
		Vector points_1d;
		sdf.topoints_1d(points_1d);
		SparseMatrix data_1d_mat;
		TensorBSplines::make_data_1d_mat(points_1d, data_1d_mat);

		SparseMatrix smooth_1d_mat; 
		TensorBSplines::make_smooth_1d_mat(smooth_1d_mat);

		SparseMatrix left_mat = data_1d_mat.transpose() * data_1d_mat + smooth_1d_mat * lambda * grid_size;
		Eigen::ConjugateGradient<SparseMatrix> solver;
		solver.compute(left_mat);

		std::cerr << "\nbegin separate direction solving:" << std::endl;

		std::cerr << "begin solving x direction:" << std::endl;
		std::cerr << "equation number = " << grid_size*grid_size << ", equation size = " << grid_size << " * " << N << std::endl;
		Vector values1(grid_size*grid_size*N);
		for (Index k = 0; k < grid_size; k++) {
			for (Index j = 0; j < grid_size; j++) {
				Index start_index = k*grid_size*grid_size + j*grid_size;
				Vector right_vec = data_1d_mat.transpose() * values.middleRows(start_index, grid_size);
				values1.middleRows(k*grid_size*N+j*N, N) = solver.solve(right_vec); 
			}
		}
		std::cerr << "finished solving x direction" << std::endl;

		std::cerr << "begin solving y direction:" << std::endl;
		std::cerr << "equation number = " << grid_size*N << ", equation size = " << grid_size << " * " << N << std::endl;
		Vector values2(grid_size*N*N);
		for (Index k = 0; k < grid_size; k++) {
			Eigen::Map<Matrix> map_values1(values1.data()+k*grid_size*N, N, grid_size);
			Matrix map_values1_traspose = map_values1.transpose();
			for (Index h = 0; h < N; h++) {
				Vector right_vec = data_1d_mat.transpose() * map_values1_traspose.col(h);
				values2.middleRows(k*N*N+h*N, N) = solver.solve(right_vec); 
			}
		}
		std::cerr << "finished solving y direction" << std::endl;

		std::cerr << "begin solving z direction:" << std::endl;
		std::cerr << "equation number = " << N*N << ", equation size = " << grid_size << " * " << N << std::endl;
		Vector values3(N*N*N);
		Eigen::Map<Matrix> map_values2(values2.data(), N*N, grid_size);
		Matrix map_values2_transpose = map_values2.transpose();
		for (Index i = 0; i < N; i++) {
			for (Index j = 0; j < N; j++) {
				Vector right_vec = data_1d_mat.transpose() * map_values2_transpose.col(i*N+j);
				values3.middleRows(i*N*N+j*N, N) = solver.solve(right_vec); 
			}
		}
		std::cerr << "finished solving z direction" << std::endl;

		tbs.controls_.swap(values3);

		std::cerr << "finished solving" << std::endl;
	}
}
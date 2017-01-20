#include <activetbsfitting.h>
#include <iostream>

namespace absr {
	template<Degree deg>
	void ActiveTBSFitting3L::fitting(ActiveTBS<deg> &atbs) {
		const PointSet &points = points_;
		const NormalSet &normals = normals_;
		const Scalar lambda = lambda_;
		const Scalar epsilon = epsilon_;
		const Size N = atbs.N_;
		const Size npts = (Size) points.rows();

		MapType &amp = atbs.amp_;
		Vector &controls = atbs.controls_;

		PointSet points_3L(npts*3, 3);
		points_3L << points, points+epsilon*normals, points-epsilon*normals;
		Vector values_3L(npts*3);
		values_3L << Vector::Zero(npts), Vector::Ones(npts), -Vector::Ones(npts);

		ActiveTBS<deg>::generate_active_map(amp, points_3L, N);

		SparseMatrix data_mat;
		ActiveTBS<deg>::make_data_mat(amp, points_3L, data_mat, N);

		SparseMatrix smooth_mat; 
		ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, N);

		SparseMatrix left_mat = data_mat.transpose() * data_mat + smooth_mat * lambda * npts * 3;
		Vector right_vec = data_mat.transpose() * values_3L;

		std::cerr << "\nbegin solving active 3L:" << std::endl;
		std::cerr << "equation number = 1, equation size = " << 3*npts << " * " << amp.size() << std::endl;

		Eigen::ConjugateGradient<SparseMatrix> cg;
		cg.compute(left_mat);
		controls = cg.solve(right_vec);

		std::cerr << "finished solving active 3L" << std::endl;
	}

	template<Degree deg>
	void ActiveTBSFittingJuttler::fitting(ActiveTBS<deg> &atbs) {
		const PointSet &points = points_;
		const NormalSet &normals = normals_;
		const Scalar lambda = lambda_;
		const Scalar kappa = kappa_;
		const Size N = atbs.N_;
		const Size npts = (Size) points.rows();

		MapType &amp = atbs.amp_;
		Vector &controls = atbs.controls_;

		ActiveTBS<deg>::generate_active_map(amp, points, N);

		SparseMatrix data_mat, du_mat, dv_mat, dw_mat;
		ActiveTBS<deg>::make_data_duvw_mat(amp, points, data_mat, du_mat, dv_mat, dw_mat, N);

		SparseMatrix smooth_mat; 
		ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, N);

		SparseMatrix left_mat = data_mat.transpose() * data_mat + 
			(du_mat.transpose() * du_mat + 
			dv_mat.transpose() * dv_mat + 
			dw_mat.transpose() * dw_mat) * kappa +
			smooth_mat * lambda * npts;
		Vector right_vec = (du_mat.transpose() * normals.col(0) +
			dv_mat.transpose() * normals.col(1) +
			dw_mat.transpose() * normals.col(2)) * kappa;

		std::cerr << "\nbegin solving active Juttler" << std::endl;
		std::cerr << "equation number = 1, equation size = " << npts << " * " << amp.size() << std::endl;

		Eigen::ConjugateGradient<SparseMatrix> cg;
		cg.compute(left_mat);
		controls = cg.solve(right_vec);

		std::cerr << "finished solving active Juttler" << std::endl;
	}

	template<Degree deg>
	void ActiveTBSFittingSDF::fitting(ActiveTBS<deg> &atbs) {
		const PointSet &points = points_;
		const SDF &sdf = sdf_;
		const Scalar lambda = lambda_;
		const Size N = atbs.N_;
		const Size npts = (Size) points.rows();

		MapType &amp = atbs.amp_;
		Vector &controls = atbs.controls_;

		const Vector &values_extra = sdf.values_;
		PointSet points_extra;
		sdf_.topoints(points_extra);
		Scalar extra_kappa = 0.01;

		ActiveTBS<deg>::generate_active_map(amp, points, N);

		SparseMatrix data_mat;
		ActiveTBS<deg>::make_data_mat(amp, points, data_mat, N);

		SparseMatrix extra_data_mat;
		ActiveTBS<deg>::make_data_mat(amp, points_extra, extra_data_mat, N);

		SparseMatrix smooth_mat; 
		ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, N);

		SparseMatrix left_mat = data_mat.transpose() * data_mat +
			extra_data_mat.transpose() * extra_data_mat * extra_kappa +
			smooth_mat * lambda * npts;

		Vector right_vec = extra_data_mat.transpose() * values_extra * extra_kappa;

		std::cerr << "\nbegin solving active SDF" << std::endl;
		std::cerr << "equation number = 1, equation size = " << npts << " * " << amp.size() << std::endl;

		Eigen::ConjugateGradient<SparseMatrix> cg;
		cg.compute(left_mat);
		controls = cg.solve(right_vec);

		std::cerr << "finished solving active SDF" << std::endl;

	}
};


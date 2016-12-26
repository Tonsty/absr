#include <activetbs.h>
#include <hierarchicaltbsfitting.h>
#include <iostream>

namespace absr {
	 template<Degree deg>
	void HierarchicalTBSFitting3L::fitting(HierarchicalTBS<deg> &htbs) {
		const PointSet &points = points_;
		const NormalSet &normals = normals_;
		const Scalar lambda = lambda_;
		const Scalar epsilon = epsilon_;
		const Size L = htbs.L_;
		const Size root_N = htbs.root_N_;
		const Size npts = (Size) points.rows();

		Vector &controls = htbs.controls_;

		PointSet points_3L(npts*3, 3);
		points_3L << points, points+epsilon*normals, points-epsilon*normals;
		Vector values_3L(npts*3);
		values_3L << Vector::Zero(npts), Vector::Ones(npts), -Vector::Ones(npts);

		std::cerr << "\nbegin solving active 3L:" << std::endl;

		htbs.amps_.resize(L);
		for (Index level = 0; level < L; level++) {
			MapType &amp = htbs.amps_[level];
			const Size current_N = (root_N-deg)*(1<<level) + deg;
			ActiveTBS<deg>::generate_active_map(amp, points_3L, current_N);

			SparseMatrix data_mat;
			ActiveTBS<deg>::make_data_mat(amp, points_3L, data_mat, current_N);

			SparseMatrix smooth_mat; 
			ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, current_N);

			SparseMatrix left_mat = data_mat.transpose() * data_mat + smooth_mat * lambda * npts * 3;
			Vector right_vec = data_mat.transpose() * values_3L;

			std::cerr << "equation number = 1, equation size = " << 3*npts << " * " << amp.size() << std::endl;

			Eigen::ConjugateGradient<SparseMatrix> cg;
			cg.compute(left_mat);
			Vector constrols_this_level = cg.solve(right_vec);
			Size current_controls_size = (Size) controls.size();
			controls.conservativeResize(current_controls_size + amp.size());
			controls.bottomRows(amp.size()) = constrols_this_level;

			values_3L -= data_mat * constrols_this_level;

			std::cerr << "level = " << level << ", rms residual error = " << std::sqrt(values_3L.squaredNorm()/(3*npts)) << std::endl;
		}

		std::cerr << "finished solving hierarchical 3L" << std::endl;
	}

	 template<Degree deg>
	void HierarchicalTBSFittingJuttler::fitting(HierarchicalTBS<deg> &htbs) {
		PointSet &points = points_;
		NormalSet &normals = normals_;
		const Scalar lambda = lambda_;
		const Scalar kappa = kappa_;
		const Size L = htbs.L_;
		const Size root_N = htbs.root_N_;
		const Size npts = (Size) points.rows();

		Vector &controls = htbs.controls_;
		Vector values = Vector::Zero(npts), nvalues = normals;

		std::cerr << "\nbegin solving hierarchical Juttler" << std::endl;

		htbs.amps_.resize(L);
		for (Index level = 0; level < L; level++) {
			MapType &amp = htbs.amps_[level]; 
			Vector &controls = htbs.controls_;
			const Size current_N = (root_N-deg)*(1<<level) + deg;
			ActiveTBS<deg>::generate_active_map(amp, points, current_N);

			SparseMatrix data_mat, du_mat, dv_mat, dw_mat;
			ActiveTBS<deg>::make_data_duvw_mat(amp, points, data_mat, du_mat, dv_mat, dw_mat, current_N);

			SparseMatrix smooth_mat; 
			ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, current_N);

			SparseMatrix left_mat = data_mat.transpose() * data_mat + 
				(du_mat.transpose() * du_mat + 
				dv_mat.transpose() * dv_mat + 
				dw_mat.transpose() * dw_mat) * kappa +
				smooth_mat * lambda * npts;
			Vector right_vec = data_mat.transpose() * values + (du_mat.transpose() * nvalues.col(0) +
				dv_mat.transpose() * nvalues.col(1) +
				dw_mat.transpose() * nvalues.col(2)) * kappa;

			std::cerr << "equation number = 1, equation size = " << npts << " * " << amp.size() << std::endl;

			Eigen::ConjugateGradient<SparseMatrix> cg;
			cg.compute(left_mat);
			Vector constrols_this_level = cg.solve(right_vec);	
			Size current_controls_size = (Size) controls.size();
			controls.conservativeResize(current_controls_size + amp.size());
			controls.bottomRows(amp.size()) = constrols_this_level;

			values -= data_mat * constrols_this_level;
			nvalues.col(0) -= du_mat * constrols_this_level;
			nvalues.col(1) -= dv_mat * constrols_this_level;
			nvalues.col(2) -= dw_mat * constrols_this_level;

			std::cerr << "level = " << level << ", rms residual error = " << std::sqrt( (values.squaredNorm()+kappa*nvalues.squaredNorm())/(3*npts) ) << std::endl;
		}

		std::cerr << "finished solving hierarchical Juttler" << std::endl;
	}
};
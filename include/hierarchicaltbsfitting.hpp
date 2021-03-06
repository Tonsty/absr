#include <activetbs.h>
#include <hierarchicaltbsfitting.h>
#include <iostream>
#include <iomanip>

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
		
		MapVecType &amps = htbs.amps_;
		Vector &controls = htbs.controls_;

		PointSet points_3L(npts*3, 3);
		points_3L << points, points+epsilon*normals, points-epsilon*normals;
		Vector values_3L(npts*3);
		values_3L << Vector::Zero(npts), Vector::Ones(npts), -Vector::Ones(npts);

		std::cerr << "\nbegin solving active 3L:" << std::endl;

		std::vector<Vector> ti0s(L), ti1s(L), ti2s(L);
		for(Index level = 0; level < L; level++) HierarchicalTBS<deg>::precompute_ti0ti1ti2(level, ti0s[level], ti1s[level], ti2s[level]);

		std::vector<Triplet> global_smooth_triplelist;
		amps.resize(L);
		for (Index level = 0; level < L; level++) {
			MapType &amp = amps[level];
			const Size current_N = (root_N-deg)*(1<<level) + deg;
			ActiveTBS<deg>::generate_active_map(amp, points_3L, current_N);

			SparseMatrix data_mat;
			ActiveTBS<deg>::make_data_mat(amp, points_3L, data_mat, current_N);

			SparseMatrix smooth_mat;
			Vector smooth_vec;
			HierarchicalTBS<deg>::make_smooth_mat_vec(amps, global_smooth_triplelist, controls, level, 
				smooth_mat, smooth_vec, ti0s, ti1s, ti2s, root_N);

			SparseMatrix left_mat = data_mat.transpose() * data_mat 
				+ smooth_mat * lambda * npts * 3;
			Vector right_vec = data_mat.transpose() * values_3L 
				+ smooth_vec * lambda * npts * 3;

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
		
		MapVecType &amps = htbs.amps_;
		Vector &controls = htbs.controls_;
		Vector values = Vector::Zero(npts);
		NormalSet nvalues = normals;

		std::cerr << "\nbegin solving hierarchical Juttler" << std::endl;

		std::vector<Vector> ti0s(L), ti1s(L), ti2s(L);
		for(Index level = 0; level < L; level++) {
			HierarchicalTBS<deg>::precompute_ti0ti1ti2(level, ti0s[level], ti1s[level], ti2s[level]);
			//std::cout << "level = " << level << std::endl;
			//for(Index h = 0; h < ti0s[level].size(); h++) 
			//	std::cout << ti0s[level](h) << "\t" << ti1s[level](h) << "\t" << ti2s[level](h) << "\t" << std::endl;
		}

		std::vector<Triplet> global_smooth_triplelist;
		amps.resize(L);
		for (Index level = 0; level < L; level++) {
			MapType &amp = amps[level]; 
			const Size current_N = ((root_N-deg)<<level) + deg;
			ActiveTBS<deg>::generate_active_map(amp, points, current_N);

			SparseMatrix data_mat, du_mat, dv_mat, dw_mat;
			ActiveTBS<deg>::make_data_duvw_mat(amp, points, data_mat, du_mat, dv_mat, dw_mat, current_N);

			SparseMatrix smooth_mat;
			Vector smooth_vec = Vector::Zero(amp.size());
			//HierarchicalTBS<deg>::make_smooth_mat_vec(amps, global_smooth_triplelist, controls, level, 
			//	smooth_mat, smooth_vec, ti0s, ti1s, ti2s, root_N);
			ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, current_N);

			//std::cout << "smooth_mat = \n" << smooth_mat.toDense() << std::endl;
			//std::cout << "smooth_vec = \n" << smooth_vec << std::endl;

			SparseMatrix left_mat = data_mat.transpose() * data_mat + 
				(du_mat.transpose() * du_mat + 
				dv_mat.transpose() * dv_mat + 
				dw_mat.transpose() * dw_mat) * kappa +
				smooth_mat * lambda * npts;
			Vector right_vec = data_mat.transpose() * values + 
				(du_mat.transpose() * nvalues.col(0) + dv_mat.transpose() * nvalues.col(1) + dw_mat.transpose() * nvalues.col(2)) * kappa + 
				smooth_vec * lambda * npts;

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

			std::cerr << "level = " << level << ", rms residual error = " << std::sqrt( (values.squaredNorm()+kappa*nvalues.squaredNorm())/npts ) << std::endl;
		}

		std::cerr << "finished solving hierarchical Juttler" << std::endl;
	}

	template<Degree deg>
	void HierarchicalTBSFittingSDF::fitting(HierarchicalTBS<deg> &htbs) {
		const PointSet &points = points_;
		const SDF &sdf = sdf_;
		const Scalar lambda = lambda_;
		const Size root_N = htbs.root_N_;
		const Size L = htbs.L_;
		const Size npts = (Size) points.rows();

		MapVecType &amps = htbs.amps_;
		Vector &controls = htbs.controls_;
		Vector values = Vector::Zero(npts);

		Vector values_extra = sdf.values_;
		PointSet points_extra;
		sdf_.topoints(points_extra);
		Scalar extra_kappa = 0.01;

		std::cerr << "\nbegin solving hierarchical SDF" << std::endl;

		std::vector<Vector> ti0s(L), ti1s(L), ti2s(L);
		for(Index level = 0; level < L; level++) {
			HierarchicalTBS<deg>::precompute_ti0ti1ti2(level, ti0s[level], ti1s[level], ti2s[level]);
			std::cerr << "level = " << level << std::endl;
			for(Index h = 0; h < ti0s[level].size(); h++) 
				std::cerr << std::setiosflags(std::ios::right) << std::setw(10) << 
				ti0s[level](h) << "\t" << ti1s[level](h) << "\t" << ti2s[level](h) << "\t" << std::endl;
		}

		std::vector<Triplet> global_smooth_triplelist;
		amps.resize(L);
		for (Index level = 0; level < L; level++) {
			MapType &amp = amps[level]; 
			const Size current_N = ((root_N-deg)<<level) + deg;
			ActiveTBS<deg>::generate_active_map(amp, points, current_N);

			SparseMatrix data_mat;
			ActiveTBS<deg>::make_data_mat(amp, points, data_mat, current_N);

			SparseMatrix extra_data_mat;
			ActiveTBS<deg>::make_data_mat(amp, points_extra, extra_data_mat, current_N);

			SparseMatrix smooth_mat;
			Vector smooth_vec = Vector::Zero(amp.size());
			HierarchicalTBS<deg>::make_smooth_mat_vec(amps, global_smooth_triplelist, controls, level, 
				smooth_mat, smooth_vec, ti0s, ti1s, ti2s, root_N);
			//ActiveTBS<deg>::make_smooth_mat(amp, smooth_mat, current_N);

			//std::cerr << "smooth_mat = \n" << smooth_mat.toDense() << std::endl;
			//std::cerr << "smooth_vec = \n" << smooth_vec << std::endl;

			SparseMatrix left_mat = data_mat.transpose() * data_mat + 
				extra_data_mat.transpose() * extra_data_mat * extra_kappa +
				smooth_mat * lambda * npts;
			Vector right_vec = data_mat.transpose() * values + 
				extra_data_mat.transpose() * values_extra * extra_kappa +
				smooth_vec * lambda * npts;

			std::cerr << "equation number = 1, equation size = " << npts << " * " << amp.size() << std::endl;

			Eigen::ConjugateGradient<SparseMatrix> cg;
			cg.compute(left_mat);
			Vector constrols_this_level = cg.solve(right_vec);	
			Size current_controls_size = (Size) controls.size();
			controls.conservativeResize(current_controls_size + amp.size());
			controls.bottomRows(amp.size()) = constrols_this_level;

			values -= data_mat * constrols_this_level;
			values_extra -= extra_data_mat * constrols_this_level;

			std::cerr << "level = " << level << ", rms residual error = " << std::sqrt(values.squaredNorm()/npts) << std::endl;
		}

		std::cerr << "finished solving hierarchical SDF" << std::endl;
	}

};
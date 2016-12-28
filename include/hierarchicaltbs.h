#ifndef HIERARCHICALTBS_H
#define HIERARCHICALTBS_H

#include <typedefs.h>

namespace absr {
	template<Degree deg>
	struct HierarchicalTBS : Function{

		HierarchicalTBS(Size L = 5, Size root_N = (deg+1)) : L_(L), root_N_(root_N) {};

		virtual void evaluate(const PointSet &points, Vector &values);
		virtual Scalar operator()(Scalar x, Scalar y, Scalar z);

		Scalar operator()(Scalar x, Scalar y, Scalar z, const Size toplevels);
		void evaluate_toplevels(const PointSet &points, const Size toplevels, Vector &values);

		static void precompute_ti0ti1ti2(const Index bottomlevel, Vector &ti0, Vector &ti1, Vector &ti2);

		static void make_smooth_triplelist(const MapVecType &amps, const Index row_level, const Index column_level,
			const Index row_start_index, const Index column_start_index, std::vector<Triplet> &smooth_triplelist, const Size root_N);

		static void make_smooth_mat_vec(const MapVecType &amps, SparseMatrix &global_smooth_mat, const Vector &controls,
			const Index currentlevel, SparseMatrix &smooth_mat, Vector &smooth_vec, 
			std::vector<Vector> &ti0s, std::vector<Vector> &ti1s, std::vector<Vector> &ti2s, const Size root_N);

		MapVecType amps_;
		Vector controls_;
		Size L_;
		Size root_N_;
	};
};

#ifndef ABSR_PREINSTANTIATE
#include <hierarchicaltbs.hpp>
#endif

#endif
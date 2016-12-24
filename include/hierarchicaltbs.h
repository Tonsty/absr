#ifndef HIERARCHICALTBS_H
#define HIERARCHICALTBS_H

#include <typedefs.h>

namespace absr {
	struct HierarchicalTBS : Function{

		HierarchicalTBS(Size L = 5, Size root_N = 4) : L_(L), root_N_(root_N) {};

		virtual void evaluate(const PointSet &points, Vector &values);
		virtual Scalar operator()(Scalar x, Scalar y, Scalar z);

		Scalar operator()(Scalar x, Scalar y, Scalar z, const Size toplevels);
		void evaluate_toplevels(const PointSet &points, const Size toplevels, Vector &values);

		MapVecType amps_;
		Vector controls_;
		Size L_;
		Size root_N_;

		static void make_data_mat(const MapVecType &amps, const PointSet &points, SparseMatrix &data_mat, const Size L, const Size root_N = 4);
		static void make_smooth_mat(const MapVecType &amps, SparseMatrix &smooth_mat, const Size L, const Size root_N = 4);
		static void make_data_duvw_mat(const MapVecType &amps, const PointSet &points, SparseMatrix &data_mat, SparseMatrix &du_mat, 
			SparseMatrix &dv_mat, SparseMatrix &dw_mat, const Size L, const Size root_N = 4);
	};
};

#endif
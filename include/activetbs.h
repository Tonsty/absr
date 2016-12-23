#ifndef ACTIVETBS_H
#define ACTIVETBS_H

#include <unordered_map>
#include <typedefs.h>

namespace absr {
	struct ActiveTBS : public Function {
		typedef std::unordered_map<Index, Index> MapType;

		ActiveTBS(Size N = 32) : N_(N) {}

		virtual void evaluate(const PointSet &points, Vector &values);
		virtual Scalar operator()(Scalar x, Scalar y, Scalar z);

		MapType amp_;
		Vector controls_;
		Size N_;

		static void generate_active_map(MapType &amp, const PointSet &points, const Size N);
		static void make_data_mat(const MapType &amp, const PointSet &points, SparseMatrix &data_mat, const Size N);
		static void make_smooth_mat(const MapType &amp, SparseMatrix &smooth_mat, const Size N);
		static void make_data_duvw_mat(const MapType &amp, const PointSet &points, SparseMatrix &data_mat, SparseMatrix &du_mat, 
			SparseMatrix &dv_mat, SparseMatrix &dw_mat, const Size N);
	}; 
};

#endif
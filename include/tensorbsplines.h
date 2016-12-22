#ifndef ABSR_TENSORBSPLINES_H
#define ABSR_TENSORBSPLINES_H

#include <common.h>
#include <sdf.h>

namespace absr {
	struct TensorBSplines {
		TensorBSplines(Size N = 32) : N_(N) {};
		void evaluate(const PointSet &points, Vector &values);
		Size N_;
		Vector controls_;
		Scalar operator()(Scalar x, Scalar y, Scalar z);

		static void make_smooth_1d_mat(SparseMatrix &smooth_1d_mat, const Size N);
		static void make_data_1d_mat(const Vector &points_1d, SparseMatrix &data_1d_mat, const Size N);
		static void make_smooth_mat(SparseMatrix &smooth_mat, const Size N);
		static void make_data_mat(const PointSet &points, SparseMatrix &data_mat, const Size N);
		static void make_data_duvw_mat(const PointSet &points, SparseMatrix &data_mat, SparseMatrix &du_mat, 
			SparseMatrix &dv_mat, SparseMatrix &dw_mat, const Size N);
	};
};

#endif

#ifndef ABSR_TENSORBSPLINES_H
#define ABSR_TENSORBSPLINES_H

#include <common.h>
#include <sdf.h>

namespace absr {
	struct TensorBSplines {
		TensorBSplines() {};
		void evaluate(const PointSet &points, Vector &values);
		Vector controls_;

		static void make_smooth_1d_mat(SparseMatrix &smooth_1d_mat);
		static void make_data_1d_mat(const Vector &points_1d, SparseMatrix &data_1d_mat);
		static void make_smooth_mat(SparseMatrix &smooth_mat);
		static void make_data_mat(const PointSet &points, SparseMatrix &data_mat);
		static void make_data_duvw_mat(const PointSet &points, 
			SparseMatrix &data_mat, 
			SparseMatrix &du_mat, 
			SparseMatrix &dv_mat, 
			SparseMatrix &dw_mat);

	};
};



#endif

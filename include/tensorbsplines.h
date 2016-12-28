#ifndef ABSR_TENSORBSPLINES_H
#define ABSR_TENSORBSPLINES_H

#include <typedefs.h>
#include <sdf.h>
#include <math_help.h>

namespace absr {
	template<Degree deg = 2, Degree der = 0>
	struct PolyFormSpline : Function1D {
		Matrix B;
		PolyFormSpline();
		Scalar operator()(Scalar x);
	};

	template<Degree deg = 2>
	struct TensorBSplines : Function {
		TensorBSplines(Size N = 32) : N_(N) {};
		virtual void evaluate(const PointSet &points, Vector &values);
		virtual Scalar operator()(Scalar x, Scalar y, Scalar z);

		Size N_;
		Vector controls_;

		static void make_smooth_1d_mat(SparseMatrix &smooth_1d_mat, const Size N);
		static void make_data_1d_mat(const Vector &points_1d, SparseMatrix &data_1d_mat, const Size N);
		static void make_smooth_mat(SparseMatrix &smooth_mat, const Size N);
		static void make_data_mat(const PointSet &points, SparseMatrix &data_mat, const Size N);
		static void make_data_duvw_mat(const PointSet &points, SparseMatrix &data_mat, SparseMatrix &du_mat, 
			SparseMatrix &dv_mat, SparseMatrix &dw_mat, const Size N);

		static void make_bs_iws(IndexWeightVec &iws, Scalar pos, const Size N);
		static void make_tbs_iws(IndexWeightVec &iws, Scalar xpos, Scalar ypos, Scalar zpos, const Size N);
		static void make_tbs_iws_duvw_ws(IndexWeightVec &iws, std::vector<Scalar> &du_ws,
			std::vector<Scalar> &dv_ws, std::vector<Scalar> &dw_ws, Scalar xpos, Scalar ypos, Scalar zpos, const Size N);

		static void precompute_pi2(Matrix &pi2);
		static void precompute_Pi2(Matrix &Pi2, const Size N);
		static void precompute_pi0pi1pi2(Matrix &pi0, Matrix &pi1, Matrix &pi2);
		static void precompute_Pi0Pi1Pi2(Matrix &Pi0, Matrix &Pi1, Matrix &Pi2, const Size N);
		
		static PolyFormSpline<deg, 0> pfs0;
		static PolyFormSpline<deg, 1> pfs1;
		static PolyFormSpline<deg, 2> pfs2;
	};
};

#ifndef ABSR_PREINSTANTIATE
#include <tensorbsplines.hpp>
#endif

#endif

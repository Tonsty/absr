#ifndef ABSR_ABSR_H
#define ABSR_ABSR_H

#include <common.h>
#include <sdf.h>
#include <tensorbsplines.h>

namespace absr {

	class TBSFitting {
	public:
		virtual void fitting(TensorBSplines &tbs)=0;
	};

	class TBSFitting3L: public TBSFitting {
	public:
		TBSFitting3L(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_(0.07), epsilon_(0.01) {}
		void set_parameters(Scalar lambda, Scalar epsilon) {
			lambda_ = lambda;
			epsilon_ = epsilon;
		}
		virtual void fitting(TensorBSplines &tbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar epsilon_;
	};

	class TBSFittingJuttler: public TBSFitting {
	public:
		TBSFittingJuttler(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_(0.08), kappa_(0.05) {}
		void set_parameters(Scalar lambda, Scalar kappa) {
			lambda_ = lambda;
			kappa_ = kappa;
		}
		virtual void fitting(TensorBSplines &tbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar kappa_;
	};

	class TBSFittingSDF: public TBSFitting {
	public:
		TBSFittingSDF(const SDF &sdf) : sdf_(sdf), lambda_(0.1), global_fitting_(false) {}
		void set_parameters(Scalar lambda, bool global_fitting = false) {
			lambda_ = lambda;
			global_fitting_ = global_fitting;
		}
		virtual void fitting(TensorBSplines &tbs);
	protected:
		SDF sdf_;
		Scalar lambda_;
		bool global_fitting_;
	};
}

#endif
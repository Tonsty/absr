#ifndef ABSR_ABSR_H
#define ABSR_ABSR_H

#include <typedefs.h>
#include <sdf.h>
#include <tensorbsplines.h>

namespace absr {

	class TBSFitting3L {
	public:
		TBSFitting3L(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_((Scalar) 0.07), epsilon_((Scalar) 0.01) {}
		void set_parameters(Scalar lambda, Scalar epsilon) {
			lambda_ = lambda;
			epsilon_ = epsilon;
		}
		template<Degree deg>
		void fitting(TensorBSplines<deg> &tbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar epsilon_;
	};

	class TBSFittingJuttler {
	public:
		TBSFittingJuttler(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_((Scalar) 0.08), kappa_((Scalar) 0.05) {}
		void set_parameters(Scalar lambda, Scalar kappa) {
			lambda_ = lambda;
			kappa_ = kappa;
		}
		template<Degree deg>
		void fitting(TensorBSplines<deg> &tbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar kappa_;
	};

	class TBSFittingSDF {
	public:
		TBSFittingSDF(const SDF &sdf) : sdf_(sdf), lambda_((Scalar) 0.1), global_fitting_(false) {}
		void set_parameters(Scalar lambda, bool global_fitting = false) {
			lambda_ = lambda;
			global_fitting_ = global_fitting;
		}
		template<Degree deg>
		void fitting(TensorBSplines<deg> &tbs);
	protected:
		SDF sdf_;
		Scalar lambda_;
		bool global_fitting_;
	};
};

#ifndef ABSR_PREINSTANTIATE
#include <tbsfitting.hpp>
#endif

#endif
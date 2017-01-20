#ifndef ACTIVETBSFITTING_H
#define ACTIVETBSFITTING_H

#include <typedefs.h>
#include <sdf.h>
#include <activetbs.h>

namespace absr {
	
	class ActiveTBSFitting3L {
	public:
		ActiveTBSFitting3L(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_((Scalar) 0.07), epsilon_((Scalar) 0.01) {}
		  void set_parameters(Scalar lambda, Scalar epsilon) {
			  lambda_ = lambda;
			  epsilon_ = epsilon;
		  }
		  template<Degree deg>
		  void fitting(ActiveTBS<deg> &atbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar epsilon_;
	};

	class ActiveTBSFittingJuttler {
	public:
		ActiveTBSFittingJuttler(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_((Scalar) 0.08), kappa_((Scalar) 0.05) {}
		  void set_parameters(Scalar lambda, Scalar kappa) {
			  lambda_ = lambda;
			  kappa_ = kappa;
		  }
		  template<Degree deg>
		  void fitting(ActiveTBS<deg> &atbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar kappa_;
	};

	class ActiveTBSFittingSDF {
	public:
		ActiveTBSFittingSDF(const PointSet &points, const SDF &sdf) 
			: points_(points), sdf_(sdf), lambda_((Scalar) 0.1) {}
		void set_parameters(Scalar lambda) {
			lambda_ = lambda;
		}
		template<Degree deg>
		void fitting(ActiveTBS<deg> &atbs);
	protected:
		PointSet points_;
		SDF sdf_;
		Scalar lambda_;
	};
};

#ifndef ABSR_PREINSTANTIATE
#include <activetbsfitting.hpp>
#endif

#endif
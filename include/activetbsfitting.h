#ifndef ACTIVETBSFITTING_H
#define ACTIVETBSFITTING_H

#include <typedefs.h>
#include <activetbs.h>

namespace absr {
	
	class ActiveTBSFitting {
	public:
		virtual void fitting(ActiveTBS &atbs)=0;
	};

	class ActiveTBSFitting3L : public ActiveTBSFitting {
	public:
		ActiveTBSFitting3L(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_(0.07), epsilon_(0.01) {}
		  void set_parameters(Scalar lambda, Scalar epsilon) {
			  lambda_ = lambda;
			  epsilon_ = epsilon;
		  }
		  virtual void fitting(ActiveTBS &tbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar epsilon_;
	};

	class ActiveTBSFittingJuttler: public ActiveTBSFitting {
	public:
		ActiveTBSFittingJuttler(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_(0.08), kappa_(0.05) {}
		  void set_parameters(Scalar lambda, Scalar kappa) {
			  lambda_ = lambda;
			  kappa_ = kappa;
		  }
		  virtual void fitting(ActiveTBS &atbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar kappa_;
	};

};

#endif
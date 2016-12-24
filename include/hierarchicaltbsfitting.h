#ifndef HIERARCHICALTBSFITTING_H
#define HIERARCHICALTBSFITTING_H

#include <typedefs.h>
#include <hierarchicaltbs.h>

namespace absr {

	class HierarchicalTBSFitting {
	public:
		virtual void fitting(HierarchicalTBS &htbs)=0;
	};

	class HierarchicalTBSFitting3L : public HierarchicalTBSFitting {
	public:
		HierarchicalTBSFitting3L(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_((Scalar) 0.07), epsilon_((Scalar) 0.01) {}
		  void set_parameters(Scalar lambda, Scalar epsilon) {
			  lambda_ = lambda;
			  epsilon_ = epsilon;
		  }
		  virtual void fitting(HierarchicalTBS &htbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar epsilon_;
	};

	class HierarchicalTBSFittingJuttler: public HierarchicalTBSFitting {
	public:
		HierarchicalTBSFittingJuttler(const PointSet &points, const NormalSet normals) : 
		  points_(points), normals_(normals), lambda_((Scalar) 0.08), kappa_((Scalar) 0.05) {}
		  void set_parameters(Scalar lambda, Scalar kappa) {
			  lambda_ = lambda;
			  kappa_ = kappa;
		  }
		  virtual void fitting(HierarchicalTBS &htbs);
	protected:
		PointSet points_;
		NormalSet normals_;
		Scalar lambda_;
		Scalar kappa_;
	};
};

#endif
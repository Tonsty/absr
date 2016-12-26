#ifndef HIERARCHICALTBS_H
#define HIERARCHICALTBS_H

#include <typedefs.h>

namespace absr {
	template<Degree deg>
	struct HierarchicalTBS : Function{

		HierarchicalTBS(Size L = 5, Size root_N = (deg+1)) : L_(L), root_N_(root_N) {};

		virtual void evaluate(const PointSet &points, Vector &values);
		virtual Scalar operator()(Scalar x, Scalar y, Scalar z);

		Scalar operator()(Scalar x, Scalar y, Scalar z, const Size toplevels);
		void evaluate_toplevels(const PointSet &points, const Size toplevels, Vector &values);

		MapVecType amps_;
		Vector controls_;
		Size L_;
		Size root_N_;
	};
};

#ifndef ABSR_PREINSTANTIATE
#include <hierarchicaltbs.hpp>
#endif

#endif
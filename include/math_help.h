#ifndef ABSR_MATH_HELP_H
#define ABSR_MATH_HELP_H

#include <typedefs.h>

namespace absr {
	struct Function1D {
		virtual Scalar operator()(Scalar x) = 0;
	};

	struct ScaleShiftFunction1D : Function1D {
		ScaleShiftFunction1D(Function1D *original_func, Scalar scale = 1.0, Scalar shift = 0.0) 
			: original_func_(original_func), scale_(scale), shift_(shift) {}
		virtual Scalar operator()(Scalar x);
		Function1D *original_func_;
		Scalar scale_;
		Scalar shift_;
	};

	struct ProductFunction1D : Function1D {
		ProductFunction1D(Function1D *original_func1, Function1D *original_func2) 
			: original_func1_(original_func1), original_func2_(original_func2) {}
		virtual Scalar operator()(Scalar x);
		Function1D *original_func1_;
		Function1D *original_func2_;
	};

	Scalar integration(Function1D *func, Degree deg, Scalar a, Scalar b);
};

#endif
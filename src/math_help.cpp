#include <iostream>
#include <math_help.h>

namespace absr {

	Scalar ScaleShiftFunction1D::operator()(Scalar x) {
		Scalar val = original_func_->operator()(scale_*x + shift_);
		return val;
	}

	Scalar ProductFunction1D::operator ()(Scalar x) {
		Scalar val = original_func1_->operator()(x) * original_func2_->operator()(x);  
		return val;
	}

	Scalar integration(Function1D *func, Degree deg, Scalar a, Scalar b) {

		static Scalar sqrt3_d3 = (Scalar) std::sqrt(3.0)/3; 
		static Scalar sqrt_3d5 = (Scalar) std::sqrt(3.0/5); 
		static Vector coe(4), x(4);
		coe << (Scalar)(18-std::sqrt(30.0))/36, (Scalar)(18+std::sqrt(30.0))/36, 
			(Scalar)(18+std::sqrt(30.0))/36, (Scalar)(18-std::sqrt(30.0))/36;
		x << (Scalar)-std::sqrt((15+2*std::sqrt(30.0))/35), (Scalar)-std::sqrt((15-2*std::sqrt(30.0))/35), 
			(Scalar)std::sqrt((15-2*std::sqrt(30.0))/35), (Scalar)std::sqrt((15+2*std::sqrt(30.0))/35);

		ScaleShiftFunction1D nfunc(func, (b-a)/2, (b+a)/2);
		Scalar val = 0;
		if(deg <= 1) { val = 2*nfunc(0.0); }
		else if(deg <= 3) { val = nfunc(-sqrt3_d3) + nfunc(sqrt3_d3); }
		else if(deg <= 5) { val = (Scalar) 5.0/9*nfunc(-sqrt_3d5) + (Scalar) 8.0/9*nfunc(0) + (Scalar) 5.0/9*nfunc(sqrt_3d5); }
		else if(deg <= 7) { for(Index h = 0; h < 4; h++) {val += coe(h)*nfunc(x(h));} }
		else { std::cerr << "cannot calculate integration of degree larger than 7" << std::endl; exit(0); }
		return val * (b-a)/2;
	}
};



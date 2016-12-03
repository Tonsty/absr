#include <iostream>

#include <unsupported/Eigen/Splines>

#include "absr.h"
using namespace absr;

int main(int argc, char** argv) {

	Vector controls = Vector::Random(10);
	Vector knots(13);
	knots << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;

	Eigen::Spline<Scalar, 1, 2> spline(knots, controls);
	std::cout << spline.degree() << std::endl;
	std::cout << spline.basisFunctions(0.5) << std::endl;
	std::cout << spline.basisFunctionDerivatives(0.5, 0) << std::endl;

	//Eigen::SplineTraits< Eigen::Spline<Scalar, 2>  >::KnotVectorType knots(1, 2);

	//Eigen::SplineTraits< Eigen::Spline<Scalar, 2> >::ControlPointVectorType::Zero(2, 1);

	return 0;
}
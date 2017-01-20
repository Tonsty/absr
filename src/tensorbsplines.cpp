#include <typedefs.h>

namespace absr {
	Vector MULT(const Vector &U, const Vector &V) {
		Vector W = Vector::Zero(U.size() + V.size() -1);
		for(Index i = 0; i < U.size(); i++) {
			for (Index j = 0; j < V.size(); j++) {
				W(i+j) += U(i)*V(j);
			}
		}
		return W;
	}

	Matrix INTofMULT(const Matrix &A, const Matrix &B) {
		Matrix C = Matrix::Zero(A.rows(), B.rows());
		for (Index i = 0; i < A.rows(); i++) {
			for (Index j = 0; j < B.rows(); j++) {
				Vector W = MULT(A.row(i), B.row(j));
				//std::cerr << "W=\n" << W << std::endl; 
				Scalar integral_sum = 0;
				for(Index k = 0; k < W.size(); k++) {
					integral_sum += W(k)/(k+1);
				}
				C(i,j) = integral_sum;
			}
		}
		return C;
	}
};

#ifdef ABSR_PREINSTANTIATE

#include <tensorbsplines.hpp>

using namespace absr;

template struct TensorBSplines<2>;
//template struct TensorBSplines<3>;

#endif


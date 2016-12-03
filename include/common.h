#ifndef ABSR_COMMON_H
#define ABSR_COMMON_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace absr {
	typedef double Scalar;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
	typedef Vector Point;
	typedef Point Normal;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> PointSet;
	typedef PointSet NormalSet;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	typedef PointSet Vertices;
	typedef int Index;
	typedef Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> Faces;

	typedef int Size;

	typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
	typedef Eigen::Triplet<Scalar> Triplet;
};

#endif
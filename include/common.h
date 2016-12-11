#ifndef ABSR_COMMON_H
#define ABSR_COMMON_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace absr {
	typedef float Scalar;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
	typedef Vector Point;
	typedef Point Normal;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> PointSet;
	typedef PointSet NormalSet;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	typedef int Index;
	typedef Index Size;

	typedef Eigen::Matrix<Index, Eigen::Dynamic, 1> Vectori;
	typedef Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> Matrixi;

	typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> Array;
	typedef Eigen::Array<Index, Eigen::Dynamic, Eigen::Dynamic> Arrayi;

	typedef PointSet Vertices;

	typedef Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> Faces;

	typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
	typedef Eigen::Triplet<Scalar> Triplet;

	typedef Matrix TransformMat;
};

#endif
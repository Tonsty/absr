#include <boundingbox.h>
#include <iostream>

using namespace absr;

BoundingBox::BoundingBox(const PointSet &points) : low_(3), high_(3) {
	Size npts = (Size) points.rows();
	if(npts == 0) return;
	low_ = points.colwise().minCoeff();
	high_ = points.colwise().maxCoeff();
}

TransformMat BoundingBox::normalize_to_unit_cube(PointSet &points) {
	BoundingBox bbox(points);

	Scalar delta = (Scalar) 0.05 * (bbox.high_-bbox.low_).norm();
	Vector delta3(3);
	delta3 << delta, delta, delta;
	Vector low = bbox.low_ - delta3, high = bbox.high_ + delta3;

	Vector shift1 = (Scalar) -0.5*(low+high);
	Scalar scale = (Scalar) 1.0/(high - low).maxCoeff();
	Vector shift2 = (Scalar) 0.5*Vector::Ones(3);

	Vector shift = scale*shift1+shift2;
	points *= scale;
	points.rowwise() += shift.transpose();

	TransformMat tf_mat = TransformMat::Identity(4, 4);
	tf_mat.diagonal() << scale, scale, scale, 1;
	tf_mat.block(0, 3, 3, 1) = shift;

	return tf_mat;
}
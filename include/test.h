#ifndef ABSR_TEST_H
#define ABSR_TEST_H

#include <typedefs.h>
#include <tensorbsplines.h>
#include <activetbs.h>
#include <hierarchicaltbs.h>

using namespace absr;

TransformMat prepare_points_normals(PointSet &points, NormalSet &normals,
	const std::string file);

void test_sdf_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar narrow_band_width, Size sdf_grid_size, TensorBSplines<2> &tbs);

void test_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, TensorBSplines<2> &tbs);

void test_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, TensorBSplines<2> &tbs);

void test_active_sdf_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar narrow_band_width, Size sdf_grid_size, ActiveTBS<2> &atbs);

void test_active_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, ActiveTBS<2> &atbs);

void test_active_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, ActiveTBS<2> &atbs);

void test_hierarchical_sdf_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar narrow_band_width, Size sdf_grid_size, HierarchicalTBS<2> &htbs);

void test_hierarchical_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, HierarchicalTBS<2> &htbs);

void test_hierarchical_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, HierarchicalTBS<2> &htbs);

#endif
#ifndef ABSR_TEST_H
#define ABSR_TEST_H

#include <typedefs.h>
#include <tensorbsplines.h>
#include <activetbs.h>
#include <hierarchicaltbs.h>

using namespace absr;

TransformMat prepare_points_normals(PointSet &points, NormalSet &normals,
	const std::string points_file, const std::string normals_file);

void test_sdf_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar narrow_band_width, Size sdf_grid_size, TensorBSplines &tbs);

void test_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, TensorBSplines &tbs);

void test_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, TensorBSplines &tbs);

void test_active_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, ActiveTBS &atbs);

void test_active_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, ActiveTBS &atbs);

void test_hierarchical_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, HierarchicalTBS &htbs);

void test_hierarchical_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, HierarchicalTBS &htbs);

void vtk_mc_display(const SDF &sdf, Scalar isovalue = 0.0, 
	const TransformMat &transmat = TransformMat::Identity(4, 4), const bool invertface = true);

#endif
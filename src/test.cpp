#include <sdf.h>
#include <fastmarching.h>
#include <tbsfitting.h>
#include <activetbsfitting.h>
#include <hierarchicaltbsfitting.h>
#include <io.h>
#include <boundingbox.h>

#include <vtk_help.h>

extern int gsave;  

using namespace absr;

TransformMat prepare_points_normals(PointSet &points, NormalSet &normals,
	const std::string file) 
{
	IO::load_points_normals(file, points, normals);
	return BoundingBox::normalize_to_unit_cube(points);
}

//void test_sdf_fitting(const PointSet &points, const NormalSet &normals, 
//	Scalar lambda, Scalar narrow_band_width, Size sdf_grid_size, TensorBSplines<2> &tbs) 
//{	
//	SDF sdf;
//	sdf.grid_size_ = sdf_grid_size;
//	sdf.voxel_length_ = (Scalar) 1.0/(sdf.grid_size_-1);
//	FastMarching fm;
//	fm.grid_size_ = sdf.grid_size_;
//	fm.voxel_length_ = sdf.voxel_length_;
//	fm.compute(points, narrow_band_width);
//	fm.tagging();
//	sdf.values_.swap(fm.values_);
//	if(gsave) savesdf(sdf, "volume.vti");
//
//	TBSFittingSDF fitter(sdf);
//	fitter.set_parameters(lambda);
//	fitter.fitting(tbs);
//}
//
//void test_3L_fitting(const PointSet &points, const NormalSet &normals, 
//	Scalar lambda, Scalar epsilon, TensorBSplines<2> &tbs) {
//		TBSFitting3L fitter(points, normals);
//		fitter.set_parameters(lambda, epsilon);
//		fitter.fitting(tbs);
//}
//
//void test_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
//	Scalar lambda, Scalar kappa, TensorBSplines<2> &tbs) {
//		TBSFittingJuttler fitter(points, normals);
//		fitter.set_parameters(lambda, kappa);
//		fitter.fitting(tbs);
//}
//
//void test_active_3L_fitting(const PointSet &points, const NormalSet &normals, 
//	Scalar lambda, Scalar epsilon, ActiveTBS<2> &atbs) {
//		ActiveTBSFitting3L fitter(points, normals);
//		fitter.set_parameters(lambda, epsilon);
//		fitter.fitting(atbs);
//}
//
//void test_active_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
//	Scalar lambda, Scalar kappa, ActiveTBS<2> &atbs) {
//		ActiveTBSFittingJuttler fitter(points, normals);
//		fitter.set_parameters(lambda, kappa);
//		fitter.fitting(atbs);
//}

void test_hierarchical_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, HierarchicalTBS<2> &htbs) {
		HierarchicalTBSFitting3L fitter(points, normals);
		fitter.set_parameters(lambda, epsilon);
		fitter.fitting(htbs);
}

void test_hierarchical_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, HierarchicalTBS<2> &htbs) {
		HierarchicalTBSFittingJuttler fitter(points, normals);
		fitter.set_parameters(lambda, kappa);
		fitter.fitting(htbs);
}
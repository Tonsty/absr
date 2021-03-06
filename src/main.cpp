#include <test.h>
#include <vtk_help.h>

#include <iostream>
#include <sstream>

#include <polygonizer.hpp>

using namespace absr;

int gsave = 0;
int guse_vtk_mc = true;

int main(int argc, char** argv) {

	std::string file = argv[1];

	std::cerr << "input file = " << file << std::endl;

	Scalar lambda = 0, aux;
	Size N = 32, mc_grid_size = 128, sdf_grid_size = 256;
	std::stringstream ss2, ss3, ss4, ss5, ss6, ss7, ss8;

	Index method = 0;
	if(argc>2) {ss2.str(argv[2]); ss2>>method;}
	if(argc>3) {ss3.str(argv[3]); ss3>>N;}
	if(argc>4) {ss4.str(argv[4]); ss4>>mc_grid_size;}
	if(argc>5) {ss5.str(argv[5]); ss5>>lambda;} 
	if(argc>6) {ss6.str(argv[6]); ss6>>aux;}
	if(argc>7) {ss7.str(argv[7]); ss7>>gsave;}
	if(argc>8) {ss8.str(argv[8]); ss8>>sdf_grid_size;}

	Scalar narrow_band_width = (Scalar) 3.5, epsilon = (Scalar) 0.01, kappa = (Scalar) 0.05;

	std::cerr << "method = " << method << std::endl;
	std::cerr << "N = " << N << std::endl;
	std::cerr << "mc_grid_size = " << mc_grid_size << std::endl;
	std::cerr << "lambda = " << lambda << std::endl;
	if(method == 0) {
		narrow_band_width = aux;
		std::cerr << "narrow_band_width = " << narrow_band_width << std::endl;
	}
	else if (method == 1) {
		epsilon = aux;
		std::cerr << "epsilon = " << epsilon << std::endl;
	}
	else if(method == 2) {
		kappa = aux;
		std::cerr << "kappa = " << kappa << std::endl;
	}
	std::cerr << "save = " << gsave << std::endl;
	if(method == 0) std::cerr << "sdf_grid_size = " << sdf_grid_size << std::endl;

	std::string yes_no;
	std::cerr << "yes/no: ";
	std::cin >> yes_no;
	if((yes_no != "y") && (yes_no != "Y") && (yes_no != "yes") && (yes_no != "Yes")) exit(0);

	PointSet points;
	NormalSet normals;
	TransformMat transmat = prepare_points_normals(points, normals, file);
	
	//TensorBSplines<2> tbs(N);
	//switch(method) {
	//case 0 : {test_sdf_fitting(points, normals, lambda, narrow_band_width, sdf_grid_size, tbs); break;}
	//case 1 : {test_3L_fitting(points, normals, lambda, epsilon, tbs); break;}
	//case 2 : {test_Juttler_fitting(points, normals, lambda, kappa, tbs); break;}
	//}
	//Function *f = &tbs;

	//ActiveTBS<2> atbs(N);
	//switch(method) {
	//case 0 : {test_active_sdf_fitting(points, normals, lambda, narrow_band_width, sdf_grid_size, atbs); break;}
	//case 1 : {test_active_3L_fitting(points, normals, lambda, epsilon, atbs); break;}
	//case 2 : {test_active_Juttler_fitting(points, normals, lambda, kappa, atbs); break;}
	//}
	//Function *f = &atbs;

	//Vector t0, t1, t2;
	//HierarchicalTBS<2>::precompute_ti0ti1ti2(1, t0, t1, t2);
	//for (Index i = 0; i < t0.size(); i++) {
	//	std::cerr << "t0(" << i << ") = " << t0(i) << ", "
	//		<< "t1(" << i << ") = " << t1(i) << ", "
	//		<< "t2(" << i << ") = " << t2(i) << std::endl;
	//}

	Size L = N;
	HierarchicalTBS<2> htbs(L, 2+16);
	switch(method) {
	case 0 : {test_hierarchical_sdf_fitting(points, normals, lambda, narrow_band_width, sdf_grid_size, htbs); break;}
	case 1 : {test_hierarchical_3L_fitting(points, normals, lambda, epsilon, htbs); break;}
	case 2 : {test_hierarchical_Juttler_fitting(points, normals, lambda, kappa, htbs); break;}
	}
	Function *f = &htbs;

	if (guse_vtk_mc) {
		SDF mc_sdf;
		mc_sdf.grid_size_ = mc_grid_size;
		mc_sdf.voxel_length_ = (Scalar) 1.0/(mc_sdf.grid_size_-1);
		PointSet mc_points;
		mc_sdf.topoints(mc_points);
		f->evaluate(mc_points, mc_sdf.values_);

		if(method == 0) {
			TransformMat swapxzmat = TransformMat::Identity(4, 4);
			swapxzmat << 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;
			vtk_mc_display(mc_sdf, 0.0, transmat.inverse() * swapxzmat, true);
		} else vtk_mc_display(mc_sdf, 0.0, transmat.inverse(), false);
	} else {
		Polygonizer polyonizer;
		srand(46236);
		Index random_index = rand() % points.rows();
		Scalar x = points(random_index, 0), y = points(random_index, 1), z = points(random_index, 2);
		if(method == 0) {
			TransformMat swapxzmat = TransformMat::Identity(4, 4);
			swapxzmat << 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;
			polyonizer.compute(f, mc_grid_size, 0.0, x, y, z, transmat.inverse() * swapxzmat, false);
		} else polyonizer.compute(f, mc_grid_size, 0.0, x, y, z, transmat.inverse(), true);
	}

	return 0;
}
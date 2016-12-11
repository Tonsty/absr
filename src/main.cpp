#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkVoxelModeller.h>
#include <vtkSphereSource.h>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>

#include <vtkCubeSource.h>
#include <vtkProperty.h>

#include <vtkXMLImageDataWriter.h>
#include <vtkPLYWriter.h>

#include <absr.h>
#include <io.h>
#include <boundingbox.h>
#include <fastmarching.h>
#include <vtkReverseSense.h>

#include <sstream>

using namespace absr;

extern Size N;
int save = 0;

void vtk_save(const SDF &sdf) {
	vtkSmartPointer<vtkImageData> volume =
		vtkSmartPointer<vtkImageData>::New();

	Size mc_grid_size = sdf.grid_size_;
	Scalar mc_voxel_length = sdf.voxel_length_;
	volume->SetDimensions(mc_grid_size, mc_grid_size, mc_grid_size);
	volume->SetSpacing(mc_voxel_length, mc_voxel_length, mc_voxel_length);
	volume->SetOrigin(0, 0, 0);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToFloat();
	volume->AllocateScalars();
	float *ptr = (float *)volume->GetScalarPointer();
	for(Index i = 0; i < mc_grid_size*mc_grid_size*mc_grid_size; i++) {
		*(ptr++) = sdf.values_(i);
	}

	vtkSmartPointer<vtkXMLImageDataWriter> writer =
		vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("volume.vti");

	writer->SetInputConnection(volume->GetProducerPort());
	writer->Write();	
}

void vtk_mc_display(const SDF &sdf, Scalar isovalue = 0.0) {

	std::cerr << "begin marchingcubes:" << std::endl;

	vtkSmartPointer<vtkImageData> volume =
		vtkSmartPointer<vtkImageData>::New();
	double isoValue = isovalue;

	vtkSmartPointer<vtkMarchingCubes> surface = 
		vtkSmartPointer<vtkMarchingCubes>::New();

	Size mc_grid_size = sdf.grid_size_;
	Scalar mc_voxel_length = sdf.voxel_length_;
	volume->SetDimensions(mc_grid_size, mc_grid_size, mc_grid_size);
	volume->SetSpacing(mc_voxel_length, mc_voxel_length, mc_voxel_length);
	volume->SetOrigin(0, 0, 0);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToFloat();
	volume->AllocateScalars();
	float *ptr = (float *)volume->GetScalarPointer();
	for(Index i = 0; i < mc_grid_size*mc_grid_size*mc_grid_size; i++) {
		*(ptr++) = sdf.values_(i);
	}

	surface->SetInput(volume);
	surface->ComputeNormalsOn();
	surface->SetValue(0, isoValue);

	vtkSmartPointer<vtkRenderer> renderer = 
		vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(.1, .2, .3);

	vtkSmartPointer<vtkPolyDataMapper> mapper = 
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(surface->GetOutputPort());
	mapper->ScalarVisibilityOff();
	vtkSmartPointer<vtkActor> actor = 
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	renderer->AddActor(actor);

	//Unit Cube Wireframe
	{
		vtkSmartPointer<vtkCubeSource> cubeSource = 
			vtkSmartPointer<vtkCubeSource>::New();
		cubeSource->SetBounds(0, 1, 0, 1, 0, 1);
		vtkSmartPointer<vtkPolyDataMapper> mapper = 
			vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(cubeSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = 
			vtkSmartPointer<vtkActor>::New();
		actor->GetProperty()->SetRepresentationToWireframe();
		actor->SetMapper(mapper);
		actor->GetProperty()->LightingOff();
		renderer->AddActor(actor);
	}

	vtkSmartPointer<vtkRenderWindow> renderWindow = 
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> interactor = 
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	renderWindow->Render();

	std::cerr << "finished marchingcubes" << std::endl;

	if (save) {
		vtkSmartPointer<vtkReverseSense> reverse = 
			vtkSmartPointer<vtkReverseSense>::New();
		reverse->SetInputConnection(surface->GetOutputPort());
		reverse->ReverseCellsOn();
		//reverse->ReverseNormalsOn();

		vtkSmartPointer<vtkPLYWriter> plyWriter =
			vtkSmartPointer<vtkPLYWriter>::New();
		plyWriter->SetFileName("mesh.ply");
		plyWriter->SetInputConnection(reverse->GetOutputPort());
		plyWriter->Write();

		vtkSmartPointer<vtkXMLImageDataWriter> writer =
			vtkSmartPointer<vtkXMLImageDataWriter>::New();
		writer->SetFileName("volume_fitting.vti");

		writer->SetInputConnection(volume->GetProducerPort());
		writer->Write();	
	}

	interactor->Start();

	float delta = 0.002;
	while(1) {
		if(isoValue < -0.025 || isoValue > 0.025) delta = -delta;
		isoValue += delta;
		surface->SetValue(0, isoValue);
		surface->Update();
		renderWindow->Render();
	}
}

void test_sdf_fitting(SDF &mc_sdf, const std::string points_file, const std::string normals_file, Scalar lambda, Size sdf_grid_size) {
	SDF sdf;
	sdf.grid_size_ = sdf_grid_size;
	sdf.voxel_length_ = 1.0/(sdf.grid_size_-1);
	
	//PointSet points;
	//sdf.topoints(points);
	//Size npts = points.rows();
	//sdf.values_.resize(npts);

	////Sphere
	//Point center(3);
	//center << 0.5, 0.5, 0.5;
	//for (Index index = 0; index < npts; index++) {
	//	Vector pc = points.row(index) - center.transpose();
	//	sdf.values_(index) = pc.norm() - 0.3;
	//}

	////Cone
	//Vector v(2);
	//v << 0.5, 0.5;
	//for (Index index = 0; index < npts; index++) {
	//	Point pt = points.row(index);
	//	sdf.values_(index) = (pt.topRows(2)-v).norm() - abs(pt(2)-0.5);
	//}

	{
		PointSet points;
		NormalSet normals;
		//IO::load_points_normals("data/duck.points", "data/duck.normals", points, normals);
		//IO::load_points_normals("data/bun.points", "data/bun.normals", points, normals);
		//IO::load_points_normals("data/mannequin.points", "data/mannequin.normals", points, normals);
		//IO::load_points_normals("data/bunny_with_holes.xyz", "", points, normals);
		//IO::load_points_normals("data/bunny.xyz", "", points, normals);
		//IO::load_points_normals("data/dragon.xyz", "", points, normals);
		//IO::load_points_normals("data/armadillo.xyz", "", points, normals);
		//IO::load_points_normals("data/blade.xyz", "", points, normals);
		//IO::load_points_normals("data/horse.xyz", "", points, normals);
		IO::load_points_normals(points_file, normals_file, points, normals);

		BoundingBox bbox(points);
		BoundingBox::normalize_to_unit_cube(points);

		FastMarching fm;
		fm.grid_size_ = sdf.grid_size_;
		fm.voxel_length_ = sdf.voxel_length_;
		fm.compute(points);
		fm.tagging();

		sdf.values_.swap(fm.values_);
	}

	if(save) vtk_save(sdf);
	//exit(0);

	//vtk_mc_display(sdf);
	//exit(0);

	ABSR absr_sdf(sdf);
	absr_sdf.abspline_fitting_sdf(lambda);

	absr_sdf.resample_sdf(mc_sdf);
}

void test_3L_fitting(SDF &mc_sdf, const std::string points_file, const std::string normals_file, Scalar lambda) {
	PointSet points;
	NormalSet normals;
	//IO::load_points_normals("duck.points", "duck.normals", points, normals);
	//IO::load_points_normals("bunny.points", "bunny.normals", points, normals);
	//IO::load_points_normals("mannequin.points", "mannequin.normals", points, normals);
	IO::load_points_normals(points_file, normals_file, points, normals);

	BoundingBox bbox(points);
	BoundingBox::normalize_to_unit_cube(points);

	ABSR absr_3L(points, normals);
	absr_3L.abspline_fitting_3L(lambda);

	absr_3L.resample_sdf(mc_sdf);
}

void test_Juttler_fitting(SDF &mc_sdf, const std::string points_file, const std::string normals_file, Scalar lambda, Scalar kappa) {
	PointSet points;
	NormalSet normals;
	//IO::load_points_normals("duck.points", "duck.normals", points, normals);
	//IO::load_points_normals("bunny.points", "bunny.normals", points, normals);
	//IO::load_points_normals("mannequin.points", "mannequin.normals", points, normals);
	IO::load_points_normals(points_file, normals_file, points, normals);
	BoundingBox bbox(points);
	BoundingBox::normalize_to_unit_cube(points);

	ABSR absr_Juttler(points, normals);
	absr_Juttler.abspline_fitting_Juttler(lambda, kappa);

	absr_Juttler.resample_sdf(mc_sdf);
}

int main(int argc, char** argv) {

	std::string points_file = argv[1];
	std::string normals_file = argv[2];

	std::cerr << "points_file = " << points_file << std::endl;
	std::cerr << "normals_file = " << normals_file << std::endl;

	Scalar lambda = 0, kappa = 0;
	Size mc_grid_size = 128, sdf_grid_size = 256;
	std::stringstream ss3, ss4, ss5, ss6, ss7, ss8, ss9;

	Index method = 0;
	if(argc>3) {ss3.str(argv[3]); ss3>>method;}
	if(argc>4) {ss4.str(argv[4]); ss4>>N;}
	if(argc>5) {ss5.str(argv[5]); ss5>>mc_grid_size;}
	if(argc>6) {ss6.str(argv[6]); ss6>>lambda;} 
	if(argc>7) {ss7.str(argv[7]); ss7>>kappa;}
	if(argc>8) {ss8.str(argv[8]); ss8>>save;}
	if(argc>9) {ss9.str(argv[9]); ss9>>sdf_grid_size;}


	std::cerr << "method = " << method << std::endl;
	std::cerr << "N = " << N << std::endl;
	std::cerr << "mc_grid_size = " << mc_grid_size << std::endl;
	std::cerr << "lambda = " << lambda << std::endl;
	if(method == 2) std::cerr << "kappa = " << kappa << std::endl;
	std::cerr << "save = " << save << std::endl;
	if(method == 0) std::cerr << "sdf_grid_size = " << sdf_grid_size << std::endl;


	std::string yes_no;
	std::cerr << "yes/no: " << std::endl;
	cin >> yes_no;
	if((yes_no != "y") && (yes_no != "Y") && (yes_no != "yes") && (yes_no != "Yes")) exit(0);

	SDF mc_sdf;
	mc_sdf.grid_size_ = mc_grid_size;
	mc_sdf.voxel_length_ = 1.0/(mc_sdf.grid_size_-1);

	switch(method) {
	case 0 : {test_sdf_fitting(mc_sdf, points_file, normals_file, lambda, sdf_grid_size); break;}
	case 1 : {test_3L_fitting(mc_sdf, points_file, normals_file, lambda); break;}
	case 2 : {test_Juttler_fitting(mc_sdf, points_file, normals_file, lambda, kappa); break;}
	}

	if(method == 2) vtk_mc_display(mc_sdf, -0.5/sdf_grid_size);
	else vtk_mc_display(mc_sdf);

	return 0;
}
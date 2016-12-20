#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkVoxelModeller.h>
#include <vtkSphereSource.h>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleTrackballCamera.h>

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

TransformMat gtransmat;

void vtk_save(const SDF &sdf) 
{
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
	memcpy(ptr, sdf.values_.data(), mc_grid_size*mc_grid_size*mc_grid_size*sizeof(float));

	vtkSmartPointer<vtkXMLImageDataWriter> writer =
		vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("volume.vti");

	writer->SetInputConnection(volume->GetProducerPort());
	writer->Write();	
}

void vtk_mc_display(const SDF &sdf, Scalar isovalue = 0.0, const TransformMat &transmat = TransformMat::Identity(4, 4)) 
{
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
	memcpy(ptr, sdf.values_.data(), mc_grid_size*mc_grid_size*mc_grid_size*sizeof(float));

	surface->SetInput(volume);
	surface->ComputeNormalsOn();
	surface->SetValue(0, isoValue);

	vtkSmartPointer<vtkTransform> transform = 
		vtkSmartPointer<vtkTransform>::New();
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> transmatd = transmat.cast<double>();
	transform->SetMatrix(transmatd.data());

	vtkSmartPointer<vtkTransformPolyDataFilter> tsurface = 
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	tsurface->SetTransform(transform);
	tsurface->SetInputConnection(surface->GetOutputPort());

	vtkSmartPointer<vtkRenderer> renderer = 
		vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(.1, .2, .3);

	vtkSmartPointer<vtkPolyDataMapper> mapper = 
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(tsurface->GetOutputPort());
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

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

	interactor->SetInteractorStyle(style);

	std::cerr << "finished marchingcubes" << std::endl;

	if (save) {
		vtkSmartPointer<vtkReverseSense> reverse = 
			vtkSmartPointer<vtkReverseSense>::New();
		reverse->SetInputConnection(tsurface->GetOutputPort());
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

	renderWindow->Render();
	interactor->Start();
}

TransformMat prepare_points_normals(PointSet &points, NormalSet &normals,
	const std::string points_file, const std::string normals_file) 
{
	IO::load_points_normals(points_file, normals_file, points, normals);
	return BoundingBox::normalize_to_unit_cube(points);
}


void test_sdf_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar narrow_band_width, Size sdf_grid_size, TensorBSplines &tbs) 
{	
	SDF sdf;
	sdf.grid_size_ = sdf_grid_size;
	sdf.voxel_length_ = 1.0/(sdf.grid_size_-1);
	FastMarching fm;
	fm.grid_size_ = sdf.grid_size_;
	fm.voxel_length_ = sdf.voxel_length_;
	fm.compute(points, narrow_band_width);
	fm.tagging();
	sdf.values_.swap(fm.values_);
	if(save) vtk_save(sdf);

	ABSRFittingSDF fitter(sdf);
	fitter.set_parameters(lambda);
	fitter.fitting(tbs);
}

void test_3L_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar epsilon, TensorBSplines &tbs) {
	ABSRFitting3L fitter(points, normals);
	fitter.set_parameters(lambda, epsilon);
	fitter.fitting(tbs);
}

void test_Juttler_fitting(const PointSet &points, const NormalSet &normals, 
	Scalar lambda, Scalar kappa, TensorBSplines &tbs) {
	ABSRFittingJuttler fitter(points, normals);
	fitter.set_parameters(lambda, kappa);
	fitter.fitting(tbs);
}

int main(int argc, char** argv) {

	std::string points_file = argv[1];
	std::string normals_file = argv[2];

	std::cerr << "points_file = " << points_file << std::endl;
	std::cerr << "normals_file = " << normals_file << std::endl;

	Scalar lambda = 0, aux;
	Size mc_grid_size = 128, sdf_grid_size = 256;
	std::stringstream ss3, ss4, ss5, ss6, ss7, ss8, ss9;

	Index method = 0;
	if(argc>3) {ss3.str(argv[3]); ss3>>method;}
	if(argc>4) {ss4.str(argv[4]); ss4>>N;}
	if(argc>5) {ss5.str(argv[5]); ss5>>mc_grid_size;}
	if(argc>6) {ss6.str(argv[6]); ss6>>lambda;} 
	if(argc>7) {ss7.str(argv[7]); ss7>>aux;}
	if(argc>8) {ss8.str(argv[8]); ss8>>save;}
	if(argc>9) {ss9.str(argv[9]); ss9>>sdf_grid_size;}

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
	std::cerr << "save = " << save << std::endl;
	if(method == 0) std::cerr << "sdf_grid_size = " << sdf_grid_size << std::endl;

	std::string yes_no;
	std::cerr << "yes/no: " << std::endl;
	cin >> yes_no;
	if((yes_no != "y") && (yes_no != "Y") && (yes_no != "yes") && (yes_no != "Yes")) exit(0);

	PointSet points;
	NormalSet normals;
	TransformMat transmat = prepare_points_normals(points, normals, points_file, normals_file);

	TensorBSplines tbs;
	switch(method) {
	case 0 : {test_sdf_fitting(points, normals, lambda, narrow_band_width, sdf_grid_size, tbs); break;}
	case 1 : {test_3L_fitting(points, normals, lambda, epsilon, tbs); break;}
	case 2 : {test_Juttler_fitting(points, normals, lambda, kappa, tbs); break;}
	}

	SDF mc_sdf;
	mc_sdf.grid_size_ = mc_grid_size;
	mc_sdf.voxel_length_ = 1.0/(mc_sdf.grid_size_-1);
	PointSet mc_points;
	mc_sdf.topoints(mc_points);
	tbs.evaluate(mc_points, mc_sdf.values_);

	if(method == 2) vtk_mc_display(mc_sdf, -0.5/sdf_grid_size);
	else vtk_mc_display(mc_sdf);

	return 0;
}
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

#include <absr.h>
#include <io.h>
#include <boundingbox.h>

using namespace absr;

void vtk_mc_display(const SDF &sdf) {
	vtkSmartPointer<vtkImageData> volume =
		vtkSmartPointer<vtkImageData>::New();
	double isoValue = 0.0;

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
	interactor->Start();

	float delta = 0.02;
	while(1) {
		if(isoValue < -0.3 || isoValue > 0.3) delta = -delta;
		isoValue += delta;
		surface->SetValue(0, isoValue);
		surface->Update();
		renderWindow->Render();
	}
}

void test_sdf_fitting(SDF &mc_sdf) {
	SDF sdf;
	sdf.grid_size_ = 32;
	sdf.voxel_length_ = 1.0/(sdf.grid_size_-1);
	PointSet points;
	sdf.topoints(points);
	
	Size npts = points.rows();
	sdf.values_.resize(npts);
	Point center(3);
	center << 0.5, 0.5, 0.5;
	for (Index index = 0; index < npts; index++) {
		Vector pc = points.row(index).transpose() - center;
		sdf.values_(index) = pc.norm() - 0.3;
	}

	//vtk_mc_display(sdf);
	ABSR absr_sdf(sdf);
	absr_sdf.abspline_fitting_sdf();

	absr_sdf.resample_sdf(mc_sdf);
}

void test_3L_fitting(SDF &mc_sdf) {
	PointSet points;
	NormalSet normals;
	//IO::load_points_normals("duck.points", "duck.normals", points, normals);
	//IO::load_points_normals("bunny.points", "bunny.normals", points, normals);
	IO::load_points_normals("mannequin.points", "mannequin.normals", points, normals);
	BoundingBox bbox(points);
	BoundingBox::normalize_to_unit_cube(points);

	ABSR absr_3L(points, normals);
	absr_3L.abspline_fitting_3L();

	absr_3L.resample_sdf(mc_sdf);
}

void test_Juttler_fitting(SDF &mc_sdf) {
	PointSet points;
	NormalSet normals;
	//IO::load_points_normals("duck.points", "duck.normals", points, normals);
	IO::load_points_normals("bunny.points", "bunny.normals", points, normals);
	//IO::load_points_normals("mannequin.points", "mannequin.normals", points, normals);
	BoundingBox bbox(points);
	BoundingBox::normalize_to_unit_cube(points);

	ABSR absr_Juttler(points, normals);
	absr_Juttler.abspline_fitting_3L();

	absr_Juttler.resample_sdf(mc_sdf);
}

int main(int argc, char** argv) {

	Vector v(9);
	v << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Eigen::Map<Matrix, Eigen::ColMajor, Eigen::InnerStride<>> m_v(v.data()+2, 3, 1, Eigen::InnerStride<>(3));
	std::cout << m_v << std::endl;
	//m_v << 0, 0, 0, 0, 0, 0, 0, 0, 0;
	//Vector v3 = Vector::Zero(3);
	//m_v.middleRows(1,1) = v3.transpose();
	//std::cout << m_v.middleRows(1, 1) << std::endl;
	//std::cout << v << std::endl;

	SDF mc_sdf;
	mc_sdf.grid_size_ = 100;
	mc_sdf.voxel_length_ = 1.0/(mc_sdf.grid_size_-1);

	test_sdf_fitting(mc_sdf);
	//test_3L_fitting(mc_sdf);
	//test_Juttler_fitting(mc_sdf);

	vtk_mc_display(mc_sdf);

	return 0;
}
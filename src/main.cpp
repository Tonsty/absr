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

using namespace absr;

void create_sdf(SDF &sdf) {
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
}

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

int main(int argc, char** argv) {

	SDF sdf;
	sdf.grid_size_ = 16;
	sdf.voxel_length_ = 1.0/(sdf.grid_size_-1);
	create_sdf(sdf);
	//vtk_mc_display(sdf);
	ABSR absr_sdf(sdf);
	absr_sdf.abspline_fitting_sdf();

	//PointSet points;
	//NormalSet normals;
	//IO::load_points_normals("duck.points", "duck.normals", points, normals);
	//ABSR absr_3L(points, normals);
	//absr_3L.abspline_fitting_3L();

	SDF mc_sdf;
	mc_sdf.grid_size_ = 64;
	mc_sdf.voxel_length_ = 1.0/(mc_sdf.grid_size_-1);
	absr_sdf.resample_sdf(mc_sdf);
	//absr_3L.resample_sdf(mc_sdf);

	vtk_mc_display(mc_sdf);

	return 0;
}
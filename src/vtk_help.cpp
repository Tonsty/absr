#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMarchingCubes.h>
#include <vtkCubeSource.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPLYWriter.h>
#include <vtkReverseSense.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <sdf.h>

using namespace absr;

extern int gsave;  

void savevolume(const vtkSmartPointer<vtkImageData> &volume, const char* filename) {
	vtkSmartPointer<vtkXMLImageDataWriter> writer =vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename);
	writer->SetInputConnection(volume->GetProducerPort());
	writer->Write();	
}

void sdftovolume(const SDF &sdf, vtkSmartPointer<vtkImageData> &volume) {
	Size grid_size = sdf.grid_size_;
	Scalar voxel_length = sdf.voxel_length_;
	volume->SetDimensions(grid_size, grid_size, grid_size);
	volume->SetSpacing(voxel_length, voxel_length, voxel_length);
	volume->SetOrigin(0, 0, 0);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToFloat();
	volume->AllocateScalars();
	float *ptr = (float *)volume->GetScalarPointer();
	memcpy(ptr, sdf.values_.data(), grid_size*grid_size*grid_size*sizeof(float));
}

void savesdf(const SDF &sdf, const char*filename) {
	vtkSmartPointer<vtkImageData> volume = vtkSmartPointer<vtkImageData>::New();
	sdftovolume(sdf, volume);
	savevolume(volume, filename);
}

vtkAlgorithmOutput* volumetosurface(vtkSmartPointer<vtkMarchingCubes> &mc, 
	vtkSmartPointer<vtkImageData> &volume, Scalar isovalue = 0.0) {
		mc->SetInput(volume);
		mc->ComputeNormalsOn();
		mc->SetValue(0, (double) isovalue);
		return mc->GetOutputPort();
}

vtkAlgorithmOutput* transformsurface(vtkSmartPointer<vtkTransformPolyDataFilter> &tfilter, vtkAlgorithmOutput*surface, 
	const TransformMat &transmat) {
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> transmatd = transmat.transpose().cast<double>();
		transform->SetMatrix(transmatd.data());
		tfilter->SetTransform(transform);
		tfilter->SetInputConnection(surface);
		return tfilter->GetOutputPort();
}

vtkAlgorithmOutput* reversesurface(vtkSmartPointer<vtkReverseSense> &reverse, vtkAlgorithmOutput*surface, 
	bool rcells, bool rnormals) {
		reverse->SetInputConnection(surface);
		reverse->ReverseCellsOff();
		reverse->ReverseNormalsOff();
		if(rcells) reverse->ReverseCellsOn();
		if(rnormals) reverse->ReverseNormalsOn();
		return reverse->GetOutputPort();
}

void drawsurface(vtkSmartPointer<vtkRenderer> &renderer, vtkAlgorithmOutput*surface) {
	vtkSmartPointer<vtkPolyDataMapper> mapper = 
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(surface);
	mapper->ScalarVisibilityOff();
	vtkSmartPointer<vtkActor> actor = 
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	renderer->AddActor(actor);
}

//Unit Cube Wireframe
void drawcube(vtkSmartPointer<vtkRenderer> &renderer, const TransformMat &transmat) {
	vtkSmartPointer<vtkCubeSource> cubeSource = 
		vtkSmartPointer<vtkCubeSource>::New();
	cubeSource->SetBounds(0, 1, 0, 1, 0, 1);

	vtkAlgorithmOutput* cube = cubeSource->GetOutputPort();

	vtkSmartPointer<vtkTransformPolyDataFilter> tfilter = 
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	cube = transformsurface(tfilter, cube, transmat);

	vtkSmartPointer<vtkPolyDataMapper> mapper = 
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(cube);
	vtkSmartPointer<vtkActor> actor = 
		vtkSmartPointer<vtkActor>::New();
	actor->GetProperty()->SetRepresentationToWireframe();
	actor->SetMapper(mapper);
	actor->GetProperty()->LightingOff();
	renderer->AddActor(actor);
}

void vtk_mc_display(const SDF &sdf, Scalar isovalue, 
	const TransformMat &transmat, const bool invertface) {

		vtkSmartPointer<vtkImageData> volume = vtkSmartPointer<vtkImageData>::New();
		std::cerr << "begin copying sdf to volume:" << std::endl;
		sdftovolume(sdf, volume);
		std::cerr << "finished copying sdf to volume:" << std::endl;
		if(gsave) {
			std::cerr << "begin saving volume_fitting.vti:" << std::endl;
			savevolume(volume, "volume_fitting.vti");
			std::cerr << "finished saving volume_fitting.vti:" << std::endl;
		}

		vtkAlgorithmOutput* surface;

		vtkSmartPointer<vtkMarchingCubes> mc = 
			vtkSmartPointer<vtkMarchingCubes>::New();
		surface= volumetosurface(mc, volume);

		vtkSmartPointer<vtkTransformPolyDataFilter> tfilter = 
			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		surface = transformsurface(tfilter, surface, transmat);
		vtkSmartPointer<vtkReverseSense> reverse = 
			vtkSmartPointer<vtkReverseSense>::New();
		if (invertface) surface = reversesurface(reverse, surface, true, false);

		vtkSmartPointer<vtkRenderer> renderer = 
			vtkSmartPointer<vtkRenderer>::New();
		renderer->SetBackground(.1, .2, .3);

		drawsurface(renderer, surface);
		drawcube(renderer, transmat);

		vtkSmartPointer<vtkRenderWindow> renderWindow = 
			vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);

		std::cerr << "begin marchingcubes:" << std::endl;
		renderWindow->Render();
		std::cerr << "finished marchingcubes" << std::endl;

		if (gsave) {
			std::cerr << "begin saving mesh.ply:" << std::endl;
			vtkSmartPointer<vtkPLYWriter> plyWriter =
				vtkSmartPointer<vtkPLYWriter>::New();
			plyWriter->SetFileName("mesh.ply");
			vtkSmartPointer<vtkReverseSense> reverse = 
				vtkSmartPointer<vtkReverseSense>::New();
			surface = reversesurface(reverse, surface, true, false);
			plyWriter->SetInputConnection(surface);
			plyWriter->Write();
			std::cerr << "finished saving mesh.ply:" << std::endl;
		}

		vtkSmartPointer<vtkRenderWindowInteractor> interactor = 
			vtkSmartPointer<vtkRenderWindowInteractor>::New();
		interactor->SetRenderWindow(renderWindow);
		vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
			vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview
		interactor->SetInteractorStyle(style);
		renderWindow->Render();
		interactor->Start();
}
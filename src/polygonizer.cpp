#include <polygonizer.h>
#include <polygonizer.hpp>
#include <tensorbsplines.h>
#include <iostream>

extern "C" int triangle2();
extern "C" char* polygonize(double (*function)(), double size, int bounds, 
	double x, double y, double z, int (*triproc)(), int mode);
extern "C" int gntris;
extern "C" VERTICES gvertices;
extern "C" TRIANGLES gtriangles;

using namespace absr;

Function *gf;
Scalar gisovalue = 0.0;

double func(double x, double y, double z) {
	if(x<0 || y<0 || z<0 || x>1 || y>1 || z>1) return 1.0;
	return gf->operator()(x, y, z) - gisovalue;
}

void Polygonizer::compute(Function *f, absr::Size grid_size, Scalar isovalue, Scalar x, Scalar y, Scalar z,
	absr::TransformMat transmat, bool invertface) {

	std::cerr << "begin polygonzation:" << std::endl;
	gf = f;
	char *err;
	fprintf(stdout, "ply\n");
	fprintf(stdout, "format ascii 1.0\n");
	fprintf(stdout, "comment polygonizer generated\n");
	if ((err = polygonize((double (*)())func, 1.0/grid_size, grid_size, x, y, z, triangle2, 0)) != NULL) {
		fprintf(stdout, "%s\n", err);
		exit(1);
	}
	std::cerr << "finished polygonzation:" << std::endl;

	fprintf(stdout, "element vertex %d\n", gvertices.count);
	fprintf(stdout, "property float x\n"); 
	fprintf(stdout, "property float y\n");
	fprintf(stdout, "property float z\n");
	fprintf(stdout, "property float nx\n"); 
	fprintf(stdout, "property float ny\n");
	fprintf(stdout, "property float nz\n");
	fprintf(stdout, "element face %d\n", gntris);
	fprintf(stdout, "property list uchar int vertex_indices\n");
	fprintf(stdout, "end_header\n");
	for (int i = 0; i < gvertices.count; i++) {
		VERTEX v;
		v = gvertices.ptr[i];
		absr::Point pt(4);
		pt << v.position.x, v.position.y, v.position.z, 1.0;
		pt = transmat * pt;
		absr::Normal nm(4);
		nm << v.normal.x, v.normal.y, v.normal.z, 0.0;
		nm = transmat * nm;
		fprintf(stdout, "%f %f %f %f %f %f\n",
			pt.x(), pt.y(),	pt.z(),
			nm.x(),	nm.y(),	nm.z());
	}
	for (int i = 0; i < gtriangles.count; i++) {
		TRIANGLE t;
		t = gtriangles.ptr[i];
		if(invertface) fprintf(stdout, "3 %d %d %d\n", t.i1, t.i3, t.i2);
		else fprintf(stdout, "3 %d %d %d\n", t.i1, t.i2, t.i3);
	}
}
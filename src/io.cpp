#include <vector>
#include <fstream>
#include <string>

#include <io.h>

using namespace absr;

void IO::load_points_normals(const std::string points_file, const std::string normals_file, PointSet &points, NormalSet &normals) {
	std::fstream file_pts(points_file, std::ios::in);
	Scalar px, py, pz;
	std::vector<Point> pvec;
	while(file_pts>>px>>py>>pz) {
		Point pt(3);
		pt << px,py,pz;
		pvec.push_back(pt);
	}

	std::fstream file_nms(normals_file, std::ios::in);
	Scalar nx, ny, nz;
	std::vector<Point> nvec;
	while(file_nms>>nx>>ny>>nz) {
		Normal nm(3);
		nm << nx,ny,nz;
		nvec.push_back(nm);
	}

	int npts = pvec.size();
	points.resize(npts, 3);
	normals.resize(npts, 3);

	for (int i = 0; i < npts; i++) {
		points.row(i) = pvec[i];
		normals.row(i) = nvec[i];
	}
}


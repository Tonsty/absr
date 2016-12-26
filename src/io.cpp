#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <typedefs.h>
#include <io.h>

using namespace absr;

void IO::load_points_normals(const std::string file, PointSet &points, NormalSet &normals) {
	std::fstream file_xyz(file, std::ios::in);

	std::vector<Point> pvec;
	std::vector<Normal> nvec;
	bool has_normal = true;
	if (file_xyz) {
		Scalar px, py, pz;
		std::string temp_line;
		while( std::getline(file_xyz, temp_line) ) {
			std::stringstream ss(temp_line);
			if (ss>>px>>py>>pz) {
				Point pt(3);
				pt << px,py,pz;
				pvec.push_back(pt);
				Scalar nx, ny, nz;
				if (has_normal) {
					if (ss>>nx>>ny>>nz) {
						Normal nm(3);
						nm << nx,ny,nz;
						nvec.push_back(nm);
					}else has_normal = false;
				}
			}
		}
	}
	file_xyz.close();

	Size npts = (Size) pvec.size(), nnms = (Size) nvec.size();
	if (npts > 0) {
		std::cerr << npts << " points, " << nnms << " normals" << std::endl;
		points.resize(npts, 3);
		normals.resize(npts, 3);
	}else {
		std::cerr << "no point loaded!!!" << std::endl;
		exit(0);
	}

	for (Index i = 0; i < npts; i++) {
		points.row(i) = pvec[i];
		if (i < nnms) {
			nvec[i].normalize();
			normals.row(i) = nvec[i];
		}
	}
}


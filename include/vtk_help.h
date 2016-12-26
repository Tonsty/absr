#ifndef ABSR_VTK_HELP_H
#define ABSR_VTK_HELP_H

#include <sdf.h>

void savesdf(const absr::SDF &sdf, const char*filename);

void vtk_mc_display(const absr::SDF &sdf, absr::Scalar isovalue = 0.0, 
	const absr::TransformMat &transmat = absr::TransformMat::Identity(4, 4), const bool invertface = true);

#endif
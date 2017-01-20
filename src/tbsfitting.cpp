#ifdef ABSR_PREINSTANTIATE

#include <tbsfitting.hpp>

using namespace absr;

template void TBSFitting3L::fitting<2>(TensorBSplines<2> &tbs);
template void TBSFittingJuttler::fitting<2>(TensorBSplines<2> &tbs);
template void TBSFittingSDF::fitting<2>(TensorBSplines<2> &tbs);

//template void TBSFitting3L::fitting<3>(TensorBSplines<3> &tbs);
//template void TBSFittingJuttler::fitting<3>(TensorBSplines<3> &tbs);
//template void TBSFittingSDF::fitting<3>(TensorBSplines<3> &tbs);

#endif
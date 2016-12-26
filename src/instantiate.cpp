#ifdef ABSR_PREINSTANTIATE
#include <typedefs.h>
#include <tensorbsplines.hpp>
#include <activetbs.hpp>
#include <hierarchicaltbs.hpp>
#include <tbsfitting.hpp>
#include <activetbsfitting.hpp>
#include <hierarchicaltbsfitting.hpp>

using namespace absr;

template struct TensorBSplines<2>;
template struct ActiveTBS<2>;
template struct HierarchicalTBS<2>;

//template struct TensorBSplines<3>;
//template struct ActiveTBS<3>;
//template struct HierarchicalTBS<3>;

//template void TBSFitting3L::fitting<2>(TensorBSplines<2> &tbs);
//template void TBSFittingJuttler::fitting<2>(TensorBSplines<2> &tbs);
//template void TBSFittingSDF::fitting<2>(TensorBSplines<2> &tbs);

//template void TBSFitting3L::fitting<3>(TensorBSplines<3> &tbs);
//template void TBSFittingJuttler::fitting<3>(TensorBSplines<3> &tbs);
//template void TBSFittingSDF::fitting<3>(TensorBSplines<3> &tbs);

//template void ActiveTBSFitting3L::fitting<2>(ActiveTBS<2> &atbs);
//template void ActiveTBSFittingJuttler::fitting<2>(ActiveTBS<2> &atbs);
//
//template void ActiveTBSFitting3L::fitting<3>(ActiveTBS<3> &atbs);
//template void ActiveTBSFittingJuttler::fitting<3>(ActiveTBS<3> &atbs);

template void HierarchicalTBSFitting3L::fitting<2>(HierarchicalTBS<2> &htbs);
template void HierarchicalTBSFittingJuttler::fitting<2>(HierarchicalTBS<2> &htbs);

//template void HierarchicalTBSFitting3L::fitting<3>(HierarchicalTBS<3> &htbs);
//template void HierarchicalTBSFittingJuttler::fitting<3>(HierarchicalTBS<3> &htbs);

#endif

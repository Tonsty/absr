#ifdef ABSR_PREINSTANTIATE

#include <hierarchicaltbsfitting.hpp>

using namespace absr;

template void HierarchicalTBSFitting3L::fitting<2>(HierarchicalTBS<2> &htbs);
template void HierarchicalTBSFittingJuttler::fitting<2>(HierarchicalTBS<2> &htbs);
template void HierarchicalTBSFittingSDF::fitting<2>(HierarchicalTBS<2> &htbs);

//template void HierarchicalTBSFitting3L::fitting<3>(HierarchicalTBS<3> &htbs);
//template void HierarchicalTBSFittingJuttler::fitting<3>(HierarchicalTBS<3> &htbs);
//template void HierarchicalTBSFitting3L::fitting<3>(HierarchicalTBS<3> &htbs);

#endif
#ifdef ABSR_PREINSTANTIATE

#include <activetbsfitting.hpp>

using namespace absr;

template void ActiveTBSFitting3L::fitting<2>(ActiveTBS<2> &atbs);
template void ActiveTBSFittingJuttler::fitting<2>(ActiveTBS<2> &atbs);
template void ActiveTBSFittingSDF::fitting(ActiveTBS<2> &atbs);

//template void ActiveTBSFitting3L::fitting<3>(ActiveTBS<3> &atbs);
//template void ActiveTBSFittingJuttler::fitting<3>(ActiveTBS<3> &atbs);
//template void ActiveTBSFittingSDF::fitting(ActiveTBS<3> &atbs);

#endif
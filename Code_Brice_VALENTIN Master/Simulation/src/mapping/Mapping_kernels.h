#ifndef MAPPING_KERNELS_H
#define MAPPING_KERNELS_H

#include <DataStructures.h>

//// EXTERNAL KERNELS ////

extern void TetraMapper3f_apply( const TVecTetra & map_i, const TVecVec4 & map_f, TVecVec3 & out, const TVecVec3 & in );

#endif

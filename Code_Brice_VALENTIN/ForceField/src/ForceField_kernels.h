#ifndef FORCEFIELD_KERNELS_H
#define FORCEFIELD_KERNELS_H

#include <DataStructures.h>

//// EXTERNAL KERNELS ////

extern void TetrahedronFEMForceField3f_initialize( TVecVec3 & positions0, TVecTetra & tetrahedra, TReal youngModulus, TReal poissonRatio, std::vector<Element> & elems);
extern void TetrahedronFEMForceField3f_addForce(TVecVec3 & positions0, std::vector<Element> & elems, TVecVec3 & f, const TVecVec3 & positions);
extern void TetrahedronFEMForceField3f_addDForce( double factor, std::vector<Element> & elems, TVecVec3 & df,const TVecVec3 & dx);

TReal invertMatrix44(TMat4x4& m, TMat4x4& out );

#endif

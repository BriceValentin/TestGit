#ifndef MASS_KERNELS_H
#define MASS_KERNELS_H

#include <DataStructures.h>

//// EXTERNAL KERNELS ////

extern void UniformMass3f_addForce( const TVec3 & gavity, TReal massDensity, TVecVec3 & f );
extern void UniformMass3f_accFromF( TReal mass, TVecVec3 & a, const TVecVec3 & f );
extern void UniformMass3f_addMDx( TReal mass, TVecVec3 & res, const TVecVec3 & dx );

#endif

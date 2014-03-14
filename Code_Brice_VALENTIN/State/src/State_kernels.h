#ifndef STATE_KERNELS_H
#define STATE_KERNELS_H

#include <DataStructures.h>

//// EXTERNAL KERNELS ////

extern void MechanicalObject3f_vPEq1( TVecVec3 & res, int index, const TVec3 val );
extern void MechanicalObject3f_vClear( TVecVec3 & res );
extern void MechanicalObject3f_vEqBF( TVecVec3 & res, const TVecVec3 & b, float f );
extern void MechanicalObject3f_vPEqBF( TVecVec3 & res, const TVecVec3 & b, float f );
extern void MechanicalObject3f_vOp( TVecVec3 & res, const TVecVec3 & a, const TVecVec3 & b, float f );
extern TReal MechanicalObject3f_vDot( const TVecVec3 & a, const TVecVec3 & b);

extern void FixedConstraint3f_projectResponseIndexed( const TVecInt & indices, TVecVec3 & dx );

#endif

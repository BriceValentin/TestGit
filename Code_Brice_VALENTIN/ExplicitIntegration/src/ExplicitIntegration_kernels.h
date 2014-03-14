#ifndef EXPLICITINTEGRATION_KERNELS_H
#define EXPLICITINTEGRATION_KERNELS_H

#include <DataStructures.h>

//// EXTERNAL KERNELS ////

extern void timeIntegrator_EulerExplicit(const SimulationParameters* params, FEMMesh* mesh, double h);

#endif

#ifndef IMPLICITINTEGRATION_KERNELS_H
#define IMPLICITINTEGRATION_KERNELS_H

#include <DataStructures.h>

//// EXTERNAL KERNELS ////

extern int simulation_cg_iter;
extern void timeIntegrator_EulerImplicit(const SimulationParameters* params, FEMMesh* mesh, double h);

#endif

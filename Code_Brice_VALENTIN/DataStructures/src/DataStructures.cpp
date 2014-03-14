#include <DataStructures.h>

//test git 3?

#include <string.h>
#include <iostream>
#include <time.h>
#include <unistd.h>

//// DATA ////
SimulationParameters simulation_params;
FEMMesh* fem_mesh = NULL;
std::vector<SurfaceMesh*> render_meshes;

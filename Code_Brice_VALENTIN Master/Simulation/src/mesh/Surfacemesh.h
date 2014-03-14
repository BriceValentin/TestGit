#ifndef SURFACEMESH_H
#define SURFACEMESH_H

#include <DataStructures.h>

void SurfaceMesh_updatePositions(SurfaceMesh * s_mesh,FEMMesh* inputMesh);
void SurfaceMesh_updateNormals(SurfaceMesh * s_mesh);
void SurfaceMesh_saveObj(SurfaceMesh * s_mesh, const std::string& filename, const std::string& mtlfilename);
void SurfaceMesh_init(SurfaceMesh * s_mesh,FEMMesh* inputMesh);
bool SurfaceMesh_load(SurfaceMesh* s_mesh,const std::string& filename);

#endif

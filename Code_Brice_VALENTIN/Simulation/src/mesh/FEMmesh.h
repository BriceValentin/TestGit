#ifndef MESH_H
#define MESH_H

#include <DataStructures.h>

void FEMMesh_removeIsolatedPoints(FEMMesh * mesh);
void FEMMesh_reorder(FEMMesh * mesh);
bool FEMMesh_isFixedParticle(FEMMesh * mesh,int index);
void FEMMesh_addFixedParticle(FEMMesh * mesh,int index);
void FEMMesh_removeFixedParticle(FEMMesh * mesh,int index);
void FEMMesh_init(FEMMesh * mesh,SimulationParameters* params);
void FEMMesh_reset(FEMMesh * mesh);
bool FEMMesh_save(FEMMesh * mesh,const std::string& filename);
bool FEMMesh_load(FEMMesh * mesh,const std::string& filename);
void FEMMesh_saveObj(FEMMesh * mesh,const std::string& filename, const std::string& mtlfilename);

#endif

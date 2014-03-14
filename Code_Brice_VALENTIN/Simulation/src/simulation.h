#ifndef SIMULATION_H
#define SIMULATION_H

#include <DataStructures.h>

bool simulation_preload();
bool simulation_load_fem_mesh();
void simulation_reorder_fem_mesh();
bool simulation_load_render_mesh();
bool simulation_init();

void simulation_animate();
void simulation_mapping();
void simulation_reset();
void simulation_save();
void simulation_load();

extern double simulation_time;
extern TVec3 simulation_bbox[2];
extern TVec3 simulation_center;
extern double simulation_size;



#endif

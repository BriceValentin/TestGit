#include "simulation.h"
#include <mesh/read_mesh_netgen.h>
#include <iostream>
#include <set>
#include <mesh/FEMmesh.h>
#include <mesh/Surfacemesh.h>
#include <octree/octree.h>
#include <State_kernels.h>
#include <ForceField_kernels.h>
#include <Mass_kernels.h>
#include <ExplicitIntegration_kernels.h>
#include <ImplicitIntegration_kernels.h>

double simulation_time = 0;
TVec3 simulation_bbox[2];
TVec3 simulation_center;
double simulation_size;
extern bool simulation_mapping_needed;

bool simulation_load_fem_mesh()
{
    FEMMesh* mesh = new FEMMesh;

    if (! simulation_params.fem_msh_filename.empty()) {
        if (!read_mesh_msh(simulation_params.fem_msh_filename.c_str(), mesh->positions, mesh->tetrahedra, mesh->triangles))
        {
            delete mesh;
            return false;
        }
        mesh->filename = simulation_params.fem_msh_filename;
    } else if (! simulation_params.fem_netget_filename.empty()) {
        if (!read_mesh_netgen(simulation_params.fem_netget_filename.c_str(), mesh->positions, mesh->tetrahedra, mesh->triangles))
        {
            delete mesh;
            return false;
        }
        mesh->filename = simulation_params.fem_netget_filename;
    } else if (! simulation_params.fem_vtk_filename.empty()) {
        if (!read_mesh_vtk(simulation_params.fem_vtk_filename.c_str(), mesh->positions, mesh->tetrahedra, mesh->triangles))
        {
            delete mesh;
            return false;
        }
        mesh->filename = simulation_params.fem_vtk_filename;
    } else {

        delete mesh;
        return false;
    }

    if (simulation_params.flip_tetra) {
        for (unsigned i=0;i<mesh->tetrahedra.size();i++) {
            int itmp = mesh->tetrahedra[i][2];
            mesh->tetrahedra[i][2] = mesh->tetrahedra[i][3];
            mesh->tetrahedra[i][3] = itmp;
        }
    }

    FEMMesh_removeIsolatedPoints(mesh);

    if (mesh && mesh->positions.size() > 0)
    {
        mesh->bbox.min = mesh->positions[0];
        mesh->bbox.max = mesh->positions[0];
        for (unsigned int i=1;i<mesh->positions.size();++i)
        {
            TVec3 p = mesh->positions[i];
            for (unsigned int c=0;c<p.size();++c)
                if (p[c] < mesh->bbox.min[c]) mesh->bbox.min[c] = p[c];
                else if (p[c] > mesh->bbox.max[c]) mesh->bbox.max[c] = p[c];
        }
    }

    fem_mesh = mesh;
    return true;
}

void simulation_reorder_fem_mesh()
{
    if (fem_mesh)
        FEMMesh_reorder(fem_mesh);
}

bool simulation_init()
{
    FEMMesh* mesh = fem_mesh;

    if (mesh)
    {
        simulation_bbox[0] = mesh->bbox.min;
        simulation_bbox[1] = mesh->bbox.max;
    }
    simulation_size = (simulation_bbox[1]-simulation_bbox[0]).norm();
    simulation_center = (simulation_bbox[0] + simulation_bbox[1]) * 0.5f;

    if (mesh)
    {
        FEMMesh_init(mesh,&simulation_params);

        TetrahedronFEMForceField3f_initialize( mesh->positions0, mesh->tetrahedra, mesh->youngModulus, mesh->poissonRatio, mesh->femElem);

        for (unsigned int i = 0; i < render_meshes.size(); ++i)
            SurfaceMesh_init(render_meshes[i],mesh);
    }

    return true;
}

void simulation_reset()
{
    FEMMesh* mesh = fem_mesh;
    FEMMesh_reset(mesh);
    mesh->simulation_mapping_needed = true;
    simulation_time = 0;
}

void simulation_save()
{
    FEMMesh* mesh = fem_mesh;
    if (simulation_time)
        FEMMesh_save(mesh,mesh->filename + ".state");
    simulation_mapping();
    std::string suffix = (simulation_time ? "-deformed" : "-initial");
    for (unsigned int i = 0; i < render_meshes.size(); ++i)
    {
        std::string filename(render_meshes[i]->filename, 0, render_meshes[i]->filename.size()-4);
        SurfaceMesh_saveObj(render_meshes[i],filename + suffix + ".obj", filename + suffix + ".mtl");
    }
    {
        std::string filename(mesh->filename, 0, mesh->filename.size()-5);
        FEMMesh_saveObj(mesh,filename + suffix + ".obj", filename + suffix + ".mtl");
    }
}

void simulation_load()
{
    FEMMesh* mesh = fem_mesh;
    FEMMesh_load(mesh,mesh->filename + ".state");
    mesh->simulation_mapping_needed = true;
    simulation_time = 1;
}

bool simulation_load_render_mesh()
{
    for (unsigned int a = 0; a < simulation_params.render_filenames.size(); ++a) {
        SurfaceMesh* s_mesh = new SurfaceMesh;
        SurfaceMesh_load(s_mesh,simulation_params.render_filenames[a]);
    }

    return true;
}



//// MAIN METHOD ////

void simulation_animate()
{
    FEMMesh* mesh = fem_mesh;
    if (!mesh) return;
    const int nSteps = simulation_params.subSteps;
    const double h  = simulation_params.timeStep / nSteps;
    for (unsigned int step = 0; step < (const unsigned) nSteps; ++step)
    {
        switch (simulation_params.odeSolver)
        {
        case ODE_EulerExplicit:
            timeIntegrator_EulerExplicit(&simulation_params, mesh, h);
            break;
        case ODE_EulerImplicit:
            timeIntegrator_EulerImplicit(&simulation_params, mesh, h);
            break;
        }
    }
    simulation_time += simulation_params.timeStep;
    mesh->simulation_mapping_needed = true;
}


void simulation_mapping()
{
    FEMMesh* mesh = fem_mesh;
    if (!mesh->simulation_mapping_needed) return;
    mesh->simulation_mapping_needed = false;

    if (!mesh) return;
    for (unsigned int i = 0; i < render_meshes.size(); ++i)
    {
        SurfaceMesh_updatePositions(render_meshes[i],mesh);
        SurfaceMesh_updateNormals(render_meshes[i]);
    }
}


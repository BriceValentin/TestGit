#include <mesh/Surfacemesh.h>
#include <simulation.h>
#include <mesh/read_mesh_netgen.h>
#include <State_kernels.h>
#include <iostream>
#include <set>
#include <mesh/write_mesh_obj.h>
#include <octree/octree.h>
#include <mesh/read_mesh_obj.h>
#include <mapping/Mapping_kernels.h>

//// INTERNAL METHODS ////

void SurfaceMesh_updatePositions(SurfaceMesh * s_mesh,FEMMesh* inputMesh)
{
    const std::vector<TVec3>& in = inputMesh->positions;
    std::vector<TVec3>& out = s_mesh->positions;
    if (s_mesh->map_f.size() != out.size() || s_mesh->map_i.size() != out.size()) return;

    TetraMapper3f_apply( s_mesh->map_i, s_mesh->map_f, out, in );
}


void SurfaceMesh_updateNormals(SurfaceMesh * s_mesh)
{
    s_mesh->normals.resize(s_mesh->positions.size());
    if (s_mesh->computeTangents)
        s_mesh->tangents.resize(s_mesh->positions.size());

    if (!s_mesh->computeTangents)
    {
        for (unsigned int i=0;i<s_mesh->normals.size();++i)
            s_mesh->normals[i].clear();

        for (unsigned int i=0;i<s_mesh->triangles.size();++i)
        {
            TVec3 n = cross(s_mesh->positions[s_mesh->triangles[i][1]]-s_mesh->positions[s_mesh->triangles[i][0]],
                             s_mesh->positions[s_mesh->triangles[i][2]]-s_mesh->positions[s_mesh->triangles[i][0]]);
            n.normalize();
            for (unsigned int j=0;j<3;++j)
                s_mesh->normals[s_mesh->triangles[i][j]] += n;
        }
        for (unsigned int i=0;i<s_mesh->normals.size();++i)
            s_mesh->normals[i].normalize();
    }
    else
    {
        for (unsigned int i=0;i<s_mesh->normals.size();++i)
        {
            s_mesh->normals[i].clear();
            s_mesh->tangents[i].clear();
        }
        for (unsigned int i=0;i<s_mesh->triangles.size();++i)
        {
            TVec3 A = s_mesh->positions[s_mesh->triangles[i][0]];
            TVec3 B = s_mesh->positions[s_mesh->triangles[i][1]];
            TVec3 C = s_mesh->positions[s_mesh->triangles[i][2]];
            B -= A;
            C -= A;
            TVec3 n = cross(B,C);
            n.normalize();
            TReal Au = s_mesh->texcoords[s_mesh->triangles[i][0]][0];
            TReal Bu = s_mesh->texcoords[s_mesh->triangles[i][1]][0];
            TReal Cu = s_mesh->texcoords[s_mesh->triangles[i][2]][0];
            Bu -= Au;
            Cu -= Au;
            TVec3 t = B * Cu - C * Bu;
            t.normalize();
            for (unsigned int j=0;j<3;++j)
            {
                s_mesh->normals[s_mesh->triangles[i][j]] += n;
                s_mesh->tangents[s_mesh->triangles[i][j]] += t;
            }
        }
        for (unsigned int i=0;i<s_mesh->normals.size();++i)
        {
            s_mesh->tangents[i] = cross(s_mesh->normals[i],s_mesh->tangents[i]);
            s_mesh->normals[i].normalize();
            s_mesh->tangents[i].normalize();
        }
    }
}

void SurfaceMesh_saveObj(SurfaceMesh * s_mesh, const std::string& filename, const std::string& mtlfilename)
{
    std::vector<int> smooth;
    smooth.push_back(1);
    std::vector<std::string> groups;
    groups.push_back("Group1");
    std::vector<std::string> mats;
    std::string mtlfile;
    if (!mtlfilename.empty())
    {
        std::string matname;
        if (!s_mesh->textureFilename.empty()) matname = "texturedMat";
        else if (s_mesh->color[0] < 0.5) matname = "darkMat";
        else matname = "defaultMat";
        mats.push_back(matname);
        std::size_t p = mtlfilename.rfind('/');
        mtlfile = std::string(mtlfilename, p == std::string::npos ? 0 : p+1);
        std::ofstream out(mtlfilename.c_str());
        if (out)
        {
            out << "newmtl " << matname << std::endl;
            out << "illum 2" << std::endl;
            out << "Ka " << s_mesh->color[0] << " " << s_mesh->color[1] << " " << s_mesh->color[2] << std::endl;
            out << "Kd " << s_mesh->color[0] << " " << s_mesh->color[1] << " " << s_mesh->color[2] << std::endl;
            out << "Ks " << 1 << " " << 1 << " " << 1 << std::endl;
            out << "Ns " << 20 << std::endl;
            if (!s_mesh->textureFilename.empty())
            {
                std::size_t p = s_mesh->textureFilename.rfind('/');
                out << "map_Kd " << s_mesh->textureFilename.substr(p == std::string::npos ? 0 : p+1) << std::endl;
            }
            out.close();
        }
    }
    write_mesh_obj(filename.c_str(), mtlfile.c_str(), s_mesh->positions, s_mesh->triangles, !s_mesh->texcoords.empty() ? &s_mesh->texcoords : NULL, &groups, &mats, &smooth);
}

void SurfaceMesh_init(SurfaceMesh * s_mesh,FEMMesh* inputMesh)
{
    // elements <-> particles table
    // FEMMesh -> SurfaceMesh mapping
    if (inputMesh)
    {
        std::cout << "Creating mapping between simulation mesh \"" << inputMesh->filename << "\" and surface mesh \"" << s_mesh->filename << "\"..." << std::endl;
        static std::string input_filename;
        static std::vector<Mat3x3d> bases;
        static std::vector<Vec3d> centers;
        static Octree<Vec3d> octree;
        const TVecTetra& tetras = inputMesh->tetrahedra;
        const std::vector<TVec3>& in = inputMesh->positions;
        const std::vector<TVec3>& out = s_mesh->positions;
        s_mesh->map_i.resize(out.size());
        s_mesh->map_f.resize(out.size());
        if (input_filename != inputMesh->filename || bases.size() != tetras.size()) // we have to recompute the octree and bases
        {
            input_filename = inputMesh->filename;
            std::vector< BBox<Vec3d> > bbox;
            bases.resize(tetras.size());
            centers.resize(tetras.size());
            bbox.resize(tetras.size());
            std::cout << "  Preparing tetrahedra" << std::endl;
            for (unsigned int t=0; t<tetras.size(); ++t)
            {
                Mat3x3d m, mt;
                m[0] = in[tetras[t][1]]-in[tetras[t][0]];
                m[1] = in[tetras[t][2]]-in[tetras[t][0]];
                m[2] = in[tetras[t][3]]-in[tetras[t][0]];
                mt.transpose(m);
                bases[t].invert(mt);
                centers[t] = (in[tetras[t][0]]+in[tetras[t][1]]+in[tetras[t][2]]+in[tetras[t][3]])*0.25;
                bbox[t].add(tetras[t].begin(), tetras[t].end(), in);
            }
            std::cout << "  Building octree" << std::endl;
            octree.init(bbox,2,2);
        }
        std::cout << "  Processing vertices" << std::endl;
        int outside = 0;
        std::vector<Octree<Vec3d>*> cells;
        for (unsigned int i=0;i<out.size();i++)
        {
            Vec3d pos = out[i];
            Vec3d coefs;
            int index = -1;
            double distance = 1e10;
            Octree<Vec3d>* cell = octree.findNear(pos);
            if (cell)
            {
                const std::vector<int>& elems = cell->elems();
                for (unsigned int e = 0; e < elems.size(); e++)
                {
                    unsigned int t = elems[e];
                    Vec3d v = bases[t] * (pos - in[tetras[t][0]]);
                    double d = std::max(std::max(-v[0],-v[1]),std::max(-v[2],v[0]+v[1]+v[2]-1));
                    if (d>0) d = (pos-centers[t]).norm2();
                    if (d<distance) { coefs = v; distance = d; index = t; }
                }
            }
            if (distance > 0)
            { // pos is outside of the fem mesh, find the nearest tetra

                // first let's find at least one tetra that is close, if not already found
                if (index >= 0) // we already have a close tetra, we need to look only for closer ones
                {
                    cells.clear();
                    octree.findAllAround(cells, pos, sqrt(distance)*1.5);
                    for (unsigned int ci = 0; ci < cells.size(); ++ci)
                    {
                        if (cells[ci] == cell) continue; // already processed this cell
                        const std::vector<int>& elems = cells[ci]->elems();
                        for (unsigned int e = 0; e < elems.size(); e++)
                        {
                            unsigned int t = elems[e];
                            double d = (pos-centers[t]).norm2();
                            if (d<distance)
                            {
                                coefs = bases[t] * (pos - in[tetras[t][0]]);
                                distance = d; index = t;
                            }
                        }
                    }
                }
                else
                {
                    // failsafe case (should not happen...), to be sure we do a brute-force search
                    for (unsigned int t = 0; t < tetras.size(); t++)
                    {
                        double d = (pos-centers[t]).norm2();
                        if (d<distance)
                        {
                            coefs = bases[t] * (pos - in[tetras[t][0]]);
                            distance = d; index = t;
                        }
                    }
                }
                if (index >= 0)
                {
                    ++outside;
                }
            }
            if (index >= 0)
            {
                //std::cout << "Surface vertex " << i << " mapped from tetra " << index << " with coefs " << coefs << std::endl;
                s_mesh->map_i[i][0] = tetras[index][0];  s_mesh->map_f[i][0] = (float)(1-coefs[0]-coefs[1]-coefs[2]);
                s_mesh->map_i[i][1] = tetras[index][1];  s_mesh->map_f[i][1] = (float)(coefs[0]);
                s_mesh->map_i[i][2] = tetras[index][2];  s_mesh->map_f[i][2] = (float)(coefs[1]);
                s_mesh->map_i[i][3] = tetras[index][3];  s_mesh->map_f[i][3] = (float)(coefs[2]);
            }
        }
        std::cout << "Mapping done: " << outside << " / " << out.size() << " vertices outside of simulation mesh" << std::endl;
    }
}


bool SurfaceMesh_load(SurfaceMesh* s_mesh,const std::string& filename)
{
    s_mesh->computeTangents = false;

    if (!read_mesh_obj(filename.c_str(), s_mesh->positions, s_mesh->triangles, (render_meshes.empty() ? &s_mesh->texcoords : NULL))) // we use the texture only for the first OBJ
    {
        delete s_mesh;
        return false;
    }

    s_mesh->filename = filename.c_str();
    SurfaceMesh_updateNormals(s_mesh);
    if (render_meshes.size()&1)
        s_mesh->color = TColor(0.15f, 0.15f, 0.15f, 1.0f);
    else
        s_mesh->color = TColor(0.8f, 0.8f, 0.8f, 1.0f);
    render_meshes.push_back(s_mesh);

    return true;
}


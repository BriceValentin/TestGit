#include <mesh/FEMmesh.h>
#include <simulation.h>
#include <mesh/read_mesh_netgen.h>
#include <State_kernels.h>
#include <iostream>
#include <set>
#include <mesh/write_mesh_obj.h>

//// INTERNAL METHODS ////

void FEMMesh_removeIsolatedPoints(FEMMesh * mesh)
{
    std::vector<int> old2newpos;
    // first count elements attached to each point
    old2newpos.resize(mesh->positions.size());
    for (unsigned int i=0;i<mesh->tetrahedra.size();++i)
        for (unsigned int j=0;j<4;++j)
            ++old2newpos[mesh->tetrahedra[i][j]];
    for (unsigned int i=0;i<mesh->triangles.size();++i)
        for (unsigned int j=0;j<3;++j)
            ++old2newpos[mesh->triangles[i][j]];
    // then computing new indices for points by skipping isolated points
    // we also copy positions to new indices
    unsigned int nbp = 0;
    for (unsigned int i=0;i<mesh->positions.size();++i)
    {
        if (old2newpos[i] > 0)
        {
            old2newpos[i] = nbp;
            if (nbp != i)
                mesh->positions[nbp] = mesh->positions[i];
            ++nbp;
        }
        else 
            old2newpos[i] = -1;
    }
    if (nbp == mesh->positions.size())
        return; // no isolated points -> nothing to do
    std::cout << "Removing " << mesh->positions.size() - nbp << " isolated points." << std::endl;

    mesh->positions.resize(nbp);
    for (unsigned int i=0;i<mesh->tetrahedra.size();++i)
        for (unsigned int j=0;j<mesh->tetrahedra[i].size();++j)
            mesh->tetrahedra[i][j] = old2newpos[mesh->tetrahedra[i][j]];
    for (unsigned int i=0;i<mesh->triangles.size();++i)
        for (unsigned int j=0;j<mesh->triangles[i].size();++j)
            mesh->triangles[i][j] = old2newpos[mesh->triangles[i][j]];
}

template<class T1, class T2>
bool SortPairFirstFn(const std::pair<T1,T2>& a, const std::pair<T1,T2>& b) { return a.first < b.first; }

void FEMMesh_reorder(FEMMesh * mesh)
{
    // simple reordering of the vertices using the largest dimension of the mesh
    int sort_coord = 0;
    if (mesh->bbox.max[1]-mesh->bbox.min[1] > mesh->bbox.max[sort_coord]-mesh->bbox.min[sort_coord])
        sort_coord = 1;
    if (mesh->bbox.max[2]-mesh->bbox.min[2] > mesh->bbox.max[sort_coord]-mesh->bbox.min[sort_coord])
        sort_coord = 2;
    std::cout << "Reordering particles based on " << (char)('X'+sort_coord) << std::endl;
    std::vector< std::pair<TReal,int> > sortp;
    sortp.resize(mesh->positions.size());
    for (unsigned int i=0;i<mesh->positions.size();++i)
        sortp[i] = std::make_pair(mesh->positions[i][sort_coord], i);
    std::sort(sortp.begin(),sortp.end(),SortPairFirstFn<TReal,int>);
    std::vector<int> old2newpos;
    old2newpos.resize(mesh->positions.size());
    for (unsigned int i=0;i<mesh->positions.size();++i)
        old2newpos[sortp[i].second] = i;
    std::vector<TVec3> newpos;
    newpos.resize(mesh->positions.size());
    for (unsigned int i=0;i<mesh->positions.size();++i)
        newpos[i] = mesh->positions[sortp[i].second];
    mesh->positions.swap(newpos);
    for (unsigned int i=0;i<mesh->tetrahedra.size();++i)
        for (unsigned int j=0;j<mesh->tetrahedra[i].size();++j)
            mesh->tetrahedra[i][j] = old2newpos[mesh->tetrahedra[i][j]];
    for (unsigned int i=0;i<mesh->triangles.size();++i)
        for (unsigned int j=0;j<mesh->triangles[i].size();++j)
            mesh->triangles[i][j] = old2newpos[mesh->triangles[i][j]];
    
    std::cout << "Reordering tetrahedra based on connected particles" << std::endl;
    std::vector< std::pair<int,int> > sortt;
    sortt.resize(mesh->tetrahedra.size());
    for (unsigned int i=0;i<mesh->tetrahedra.size();++i)
        sortt[i] = std::make_pair(std::min(std::min(mesh->tetrahedra[i][0],mesh->tetrahedra[i][1]),std::min(mesh->tetrahedra[i][2],mesh->tetrahedra[i][3])), i);
    std::sort(sortt.begin(),sortt.end(),SortPairFirstFn<int,int>);
    TVecTetra newtetra;
    newtetra.resize(mesh->tetrahedra.size());
    for (unsigned int i=0;i<mesh->tetrahedra.size();++i)
        newtetra[i] = mesh->tetrahedra[sortt[i].second];
    mesh->tetrahedra.swap(newtetra);
    std::cout << "Mesh reordering done" << std::endl;
}

bool FEMMesh_isFixedParticle(FEMMesh * mesh,int index)
{
    for (int i=0;i<mesh->nbFixedParticles;i++) {
        if (mesh->fixedParticles[i] == index) return true;
    }
    return false;
}

void FEMMesh_addFixedParticle(FEMMesh * mesh,int index)
{
    if (FEMMesh_isFixedParticle(mesh,index)) return;

    // for standard kernels we use an indices vector
    if ((unsigned) mesh->nbFixedParticles == mesh->fixedParticles.size()) {
        mesh->fixedParticles.push_back(index);
    } else {
        mesh->fixedParticles[mesh->nbFixedParticles] = index;
    }

    mesh->nbFixedParticles++;
    mesh->positions[index] = mesh->positions0[index];
    if (!mesh->velocity.empty()) mesh->velocity[index].clear();
}

void FEMMesh_removeFixedParticle(FEMMesh * mesh,int index)
{
    // for standard kernels we use an indices vector
    int i;
    for (i=0;i<mesh->nbFixedParticles;++i)
        if (mesh->fixedParticles[i] == index) break;

    if (i < mesh->nbFixedParticles)
    {
        //swap the last fixed id with this one
        mesh->fixedParticles[i] = mesh->fixedParticles[mesh->nbFixedParticles-1];
    }

    --mesh->nbFixedParticles;
}

void FEMMesh_init(FEMMesh * mesh,SimulationParameters* params)
{
    if (mesh->positions0.size() != mesh->positions.size())
    {
        mesh->positions0 = mesh->positions;
        mesh->velocity.resize(mesh->positions.size());
    }


    mesh->nbFixedParticles = 0;
    mesh->externalForce.index = -1;

    // Fixed
    mesh->fixedParticles.clear();
    mesh->nbFixedParticles = 0;
    for (unsigned i=0;i<params->fixedBBox.size();i++) {
        for (unsigned p=0;p<mesh->positions.size();p++) {
            if (mesh->positions[p][0] >= params->fixedBBox[i].min[0] && mesh->positions[p][0] <= params->fixedBBox[i].max[0] &&
                mesh->positions[p][1] >= params->fixedBBox[i].min[1] && mesh->positions[p][1] <= params->fixedBBox[i].max[1] &&
                mesh->positions[p][2] >= params->fixedBBox[i].min[2] && mesh->positions[p][2] <= params->fixedBBox[i].max[2]) {
                FEMMesh_addFixedParticle(mesh,p);
            }
        }
    }

    // FEM
    const int nbe = mesh->tetrahedra.size();
    const int nbBe = nbe;
    mesh->femElem.resize(nbBe);
    mesh->youngModulus = params->youngModulus;
    mesh->poissonRatio = params->poissonRatio;
    mesh->massDensity = params->massDensity;
}

void FEMMesh_reset(FEMMesh * mesh)
{
    mesh->positions = mesh->positions0;
    mesh->velocity.clear();
    mesh->velocity.resize(mesh->positions.size());
}

bool FEMMesh_save(FEMMesh * mesh,const std::string& filename)
{
    std::ofstream out(filename.c_str());
    if (!out)
    {
        std::cerr << "Cannot write to file " << filename << std::endl;
        return false;
    }
    out << mesh->positions.size() << " " << 6 << std::endl;
    for (unsigned int i=0;i<mesh->positions.size();++i)
    {
        out << mesh->positions[i] << " " << mesh->velocity[i] << "\n";
    }
    out.flush();
    out.close();
    return true;
}

bool FEMMesh_load(FEMMesh * mesh,const std::string& filename)
{
    std::ifstream in(filename.c_str());
    if (!in)
    {
        std::cerr << "Cannot open file " << filename << std::endl;
        return false;
    }
    int nbp = 0, nbc = 0;
    in >> nbp >> nbc;
    if ((unsigned) nbp != mesh->positions.size())
    {
        std::cerr << "ERROR: file " << mesh->filename << " contains " << nbp << " vertices while the mesh contains " << mesh->positions.size() << std::endl;
        return false;
    }
    if (nbc != 6)
    {
        std::cerr << "ERROR: file " << mesh->filename << " contains " << nbc << " values instead of 6" << std::endl;
        return false;
    }
    for (unsigned int i=0;i<mesh->positions.size();++i)
    {
        in >> mesh->positions[i] >> mesh->velocity[i];
    }
    in.close();

    return true;
}

void FEMMesh_saveObj(FEMMesh * mesh,const std::string& filename, const std::string& mtlfilename)
{
    std::vector<int> smooth;
    smooth.push_back(-1);
    std::vector<std::string> groups;
    groups.push_back("Tetras");
    std::vector<std::string> mats;
    std::string mtlfile;
    std::string matname[4];
    if (!mtlfilename.empty())
    {
        matname[0] = "fD";
        matname[1] = "fA";
        matname[2] = "fB";
        matname[3] = "fC";
        TColor colors[4];
        colors[0] = TColor(0,0,1,1);
        colors[1] = TColor(0,0.5f,1,1);
        colors[2] = TColor(0,1,1,1);
        colors[3] = TColor(0.5f,1,1,1);
        std::size_t p = mtlfilename.rfind('/');
        mtlfile = std::string(mtlfilename, p == std::string::npos ? 0 : p+1);
        std::ofstream out(mtlfilename.c_str());
        if (out)
        {
            for (int m=0;m<4;++m)
            {
                out << "newmtl " << matname[m] << std::endl;
                TColor color = colors[m];
                out << "illum 2" << std::endl;
                out << "Ka " << color[0] << " " << color[1] << " " << color[2] << std::endl;
                out << "Kd " << color[0] << " " << color[1] << " " << color[2] << std::endl;
                out << "Ks " << 1 << " " << 1 << " " << 1 << std::endl;
                out << "Ns " << 20 << std::endl;
            }
            out.close();
        }
    }
    std::vector<TVec3> points; points.reserve(mesh->tetrahedra.size()*4);
    TVecTriangle triangles; triangles.reserve(mesh->tetrahedra.size()*4);
    mats.reserve(mesh->tetrahedra.size()*4);
    for (unsigned int i=0;i<mesh->tetrahedra.size();++i)
    {
        TTetra t = mesh->tetrahedra[i];

        TVec3 a = mesh->positions[t[0]];
        TVec3 b = mesh->positions[t[1]];
        TVec3 c = mesh->positions[t[2]];
        TVec3 d = mesh->positions[t[3]];
        TVec3 center = (a+b+c+d)*(TReal)0.125;
        a = (a+center)*(TReal)0.666667;
        b = (b+center)*(TReal)0.666667;
        c = (c+center)*(TReal)0.666667;
        d = (d+center)*(TReal)0.666667;
        unsigned int pi = points.size();
        points.push_back(a);
        points.push_back(b);
        points.push_back(c);
        points.push_back(d);

        //glColor4f(0,0,1,1);
        //glVertex3fv(a.ptr()); glVertex3fv(b.ptr()); glVertex3fv(c.ptr());
        mats.push_back(matname[0]);
        triangles.push_back(TTriangle(pi+0, pi+1, pi+2));

        //glColor4f(0,0.5f,1,1);
        //glVertex3fv(b.ptr()); glVertex3fv(c.ptr()); glVertex3fv(d.ptr());
        mats.push_back(matname[1]);
        triangles.push_back(TTriangle(pi+1, pi+2, pi+3));

        //glColor4f(0,1,1,1);
        //glVertex3fv(c.ptr()); glVertex3fv(d.ptr()); glVertex3fv(a.ptr());
        mats.push_back(matname[2]);
        triangles.push_back(TTriangle(pi+2, pi+3, pi+0));

        //glColor4f(0.5f,1,1,1);
        //glVertex3fv(d.ptr()); glVertex3fv(a.ptr()); glVertex3fv(b.ptr());
        mats.push_back(matname[3]);
        triangles.push_back(TTriangle(pi+3, pi+0, pi+1));
    }

    write_mesh_obj(filename.c_str(), mtlfile.c_str(), points, triangles, (const TVecTexCoord*)NULL, &groups, &mats, &smooth);
}



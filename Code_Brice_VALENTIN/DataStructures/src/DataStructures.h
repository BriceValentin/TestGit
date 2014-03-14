#ifndef PARAMS_H
#define PARAMS_H

#include <Sofa/Vec.h>
#include <Sofa/Mat.h>

typedef float TReal;
typedef sofa::defaulttype::Vec<3,TReal> TVec3;
typedef sofa::defaulttype::Vec<4,TReal> TVec4;
typedef sofa::defaulttype::Mat<3,3,TReal> TMat3x3;
typedef sofa::defaulttype::Mat<4,4,TReal> TMat4x4;
typedef sofa::defaulttype::Mat<6,6,TReal> TMat6x6;
typedef sofa::defaulttype::Mat<12,6,TReal> TMat12x6;
typedef sofa::defaulttype::Mat<6,12,TReal> TMat6x12;
typedef sofa::defaulttype::Mat<12,12,TReal> TMat12x12;
typedef std::vector<TVec3> TVecVec3;
typedef std::vector<TVec4> TVecVec4;

typedef sofa::helper::fixed_array<unsigned int, 3> TTriangle;
typedef sofa::helper::fixed_array<unsigned int, 4> TTetra;
typedef std::vector<TTriangle> TVecTriangle;
typedef std::vector<TTetra> TVecTetra;

typedef sofa::defaulttype::Vec<2,float> TTexCoord;
typedef sofa::defaulttype::Vec<4,float> TColor;
typedef std::vector<TTexCoord> TVecTexCoord;

typedef sofa::defaulttype::Vec<3,float> Vec3f;
typedef sofa::defaulttype::Vec<3,double> Vec3d;
typedef sofa::defaulttype::Vec<4,float> Vec4f;
typedef sofa::defaulttype::Vec<4,double> Vec4d;
typedef sofa::defaulttype::Vec<4,int> Vec4i;
typedef sofa::defaulttype::Mat<3,3,float> TMat3x3f;
typedef sofa::defaulttype::Mat<3,3,double> Mat3x3d;

typedef std::vector<int> TVecInt;


enum TimeIntegration
{
    ODE_EulerExplicit = 0,
    ODE_EulerImplicit = 1
};

struct BoundingBox
{
    TVec3 min;
    TVec3 max;
};

struct SimulationParameters
{
    // Time integration
    double timeStep;
    int subSteps;
    TimeIntegration odeSolver;
    double rayleighMass;
    double rayleighStiffness;
    // CG Solver
    int maxIter;
    double tolerance;
    // Material properties
    double massDensity;
    // External forces
    TVec3 gravity;
    TReal youngModulus;
    TReal poissonRatio;

    std::string fem_msh_filename;
    std::string fem_netget_filename;
    std::string fem_vtk_filename;
    std::vector<std::string> render_filenames;
    std::string colorMap;
    std::string normalMap;
    std::string shaderVertex;
    std::string shaderFragment;

    std::string device_name;

    bool reorder;
    bool flip_tetra;
    std::vector<BoundingBox> fixedBBox;
};

struct Element
{   
    int vIx[4]; /// index of the four nodes
    TMat12x12 Ke; /// local stiffness matrix
};

struct ExternalForce
{
    int index;
    TVec3 value;
};

struct FEMMesh
{
    std::string filename;
    TVecTetra tetrahedra;
    TVecTriangle triangles;
    BoundingBox bbox;
    TVecVec3 positions;
    TVecVec3 velocity;
    TVecVec3 positions0; // rest positions

    TReal youngModulus;
    TReal poissonRatio;
    double massDensity;

    // Description of external forces
    // In a real application, this could be extended to a different force for each particle
    ExternalForce externalForce;

    // Description of constraints
    int nbFixedParticles;
    std::vector<int> fixedParticles;

    // Internal data and methods for simulation
    std::vector<TVec3> f; // force vector when using Euler explicit
    std::vector<TVec3> a,b; // solution and right-hand term when calling CG solver
    std::vector<TVec3> r,d,q; // temporary vectors used by CG solver

    // TetrahedronFEMForceField
    std::vector<Element> femElem;

    bool simulation_mapping_needed;
};

struct SurfaceMesh
{
    std::string filename;
    std::vector<TVec3> positions;
    TVecTriangle triangles;
    std::vector<TVec3> normals;
    TVecTexCoord texcoords;
    std::vector<TVec3> tangents;

    bool computeTangents;

    std::string textureFilename;
    TColor color;

    // Internal data to map FEM mesh deformation to this surface mesh
    std::vector<TTetra> map_i;
    std::vector<TVec4> map_f;
};


//External Data
extern SimulationParameters simulation_params;
extern FEMMesh* fem_mesh;
extern std::vector<SurfaceMesh*> render_meshes;

#endif

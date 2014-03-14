#ifndef READ_MESH_NETGEN_H
#define READ_MESH_NETGEN_H
#include <iostream>
#include <fstream>
#include <fstream>
#include<string>
#include<iostream>
#include <cstring>

template<class VecCoord, class VecTetra, class VecTri>
bool read_mesh_netgen(const char* filename, VecCoord& points, VecTetra& tetrahedra, VecTri& triangles)
{
    std::ifstream in(filename);
    if (!in)
    {
        std::cerr << "Cannot open file " << filename << std::endl;
        return false;
    }
    std::cout << "Reading file " << filename << std::endl;
    unsigned int nbp = 0;
    in >> nbp;
    points.resize(nbp);
    std::cout << nbp << " points" << std::endl;
    for (unsigned int i=0;i<nbp;++i)
    {
        in >> points[i][0] >> points[i][1] >> points[i][2];
    }

    unsigned int nbe = 0;
    in >> nbe;
    tetrahedra.resize(nbe);
    std::cout << nbe << " tetrahedra" << std::endl;
    for (unsigned int i=0;i<nbe;++i)
    {
        int domain;
        int a,b,c,d;
        in >> domain >> a >> b >> c >> d;
        tetrahedra[i][0] = a-1;
        tetrahedra[i][1] = b-1;
        tetrahedra[i][2] = c-1;
        tetrahedra[i][3] = d-1;
        for (unsigned int j=0;j<4;++j)
            if ((unsigned)tetrahedra[i][j] >= nbp)
            {
                std::cerr << "ERROR: invalid index " << tetrahedra[i][j] << " in tetrahedron " << i << std::endl;
                tetrahedra[i][j] = 0;
            }
    }

    unsigned int nbt = 0;
    in >> nbt;
    triangles.resize(nbt);
    std::cout << nbt << " triangles" << std::endl;
    for (unsigned int i=0;i<nbt;++i)
    {
        int boundary;
        int a,b,c;
        in >> boundary >> a >> b >> c;
        triangles[i][0] = a-1;
        triangles[i][1] = b-1;
        triangles[i][2] = c-1;
        for (unsigned int j=0;j<3;++j)
            if ((unsigned)triangles[i][j] >= nbp)
            {
                std::cerr << "ERROR: invalid index " << triangles[i][j] << " in triangle " << i << std::endl;
                triangles[i][j] = 0;
            }
    }

    in.close();

    return true;
}

template<class VecCoord>
bool read_msh_node(std::ifstream & in,VecCoord& points) {
    unsigned int nbp = 0;
    in >> nbp;
    points.resize(nbp);
    std::cout << nbp << " points" << std::endl;
    int idx;
    for (unsigned int i=0;i<nbp;++i)
    {
        in >> idx >> points[i][0] >> points[i][1] >> points[i][2];
    }
    std::string end;
    in >> end;
    return (end == "$ENDNOD");
}

template<class VecCoord, class VecTetra, class VecTri>
bool read_msh_element(std::ifstream & in,VecCoord& points, VecTetra& tetrahedra, VecTri& /*triangles*/) {
    unsigned int nbe = 0;
    in >> nbe;

    tetrahedra.resize(nbe);
    std::cout << nbe << " tetrahedra" << std::endl;


    for (unsigned int i=0;i<nbe;++i)
    {
        int idx,t1,t2,t3,t4;
        int a,b,c,d;
        in >> idx >> t1 >> t2 >> t3 >> t4 >> a >> b >> c >> d;

        if (t1==4 && t2==1 && t3==1 && t4==4) {
            tetrahedra[i][0] = a-1;
            tetrahedra[i][1] = b-1;
            tetrahedra[i][2] = c-1;
            tetrahedra[i][3] = d-1;
            for (unsigned int j=0;j<4;++j)
                if ((unsigned)tetrahedra[i][j] >= points.size())
                {
                    std::cerr << "ERROR: invalid index " << tetrahedra[i][j] << " in tetrahedron " << i << std::endl;
                    tetrahedra[i][j] = 0;
                }
        } else {
            //    unsigned int nbt = 0;
            //    in >> nbt;
            //    triangles.resize(nbt);
            //    std::cout << nbt << " triangles" << std::endl;
            //    for (unsigned int i=0;i<nbt;++i)
            //    {
            //        int boundary;
            //        int a,b,c;
            //        in >> boundary >> a >> b >> c;
            //        triangles[i][0] = a-1;
            //        triangles[i][1] = b-1;
            //        triangles[i][2] = c-1;
            //        for (unsigned int j=0;j<3;++j)
            //            if ((unsigned)triangles[i][j] >= nbp)
            //            {
            //                std::cerr << "ERROR: invalid index " << triangles[i][j] << " in triangle " << i << std::endl;
            //                triangles[i][j] = 0;
            //            }
            //    }

            std::cerr << "Cannot read element type " << t1 << " " << t2 << " "<< t3 << " " << t4 << std::endl;
            continue;
        }
    }


    std::string end;
    in >> end;
    return (end=="$ENDELM");
}

template<class VecCoord, class VecTetra, class VecTri>
bool read_mesh_msh(const char* filename, VecCoord& points, VecTetra& tetrahedra, VecTri& triangles)
{
    std::ifstream in(filename);
    if (!in)
    {
        std::cerr << "Cannot open file " << filename << std::endl;
        return false;
    }
    std::cout << "Reading file " << filename << std::endl;

    std::string type;
    while (in >> type) {
        if (type=="$NOD") {
            if (! read_msh_node(in,points)) return false;
        } else if (type=="$ELM") {
            if (! read_msh_element(in,points,tetrahedra,triangles)) return false;
        }
    }

    in.close();

    return true;
}


template<class VecCoord>
bool read_vtk_node(std::ifstream & in,VecCoord& points) {
    unsigned int nbp = 0;
    in >> nbp;
    points.resize(nbp);
    std::cout << nbp << " points" << std::endl;
    std::string type;
    in >> type;

    for (unsigned int i=0;i<nbp;++i)
    {
        in >> points[i][0] >> points[i][1] >> points[i][2];
    }

    return true;
}

template<class VecCoord, class VecTetra, class VecTri>
bool read_vtk_element(std::ifstream & in,VecCoord& points, VecTetra& tetrahedra, VecTri& /*triangles*/) {
    unsigned int nbe = 0;
    in >> nbe;

    tetrahedra.resize(nbe);
    std::cout << nbe << " tetrahedra" << std::endl;
    std::string num;
    in >> num;

    for (unsigned int i=0;i<nbe;++i)
    {
        int type;
        int a,b,c,d;
        in >> type >> a >> b >> c >> d;

        if (type == 4) {
            tetrahedra[i][0] = a;
            tetrahedra[i][1] = b;
            tetrahedra[i][2] = c;
            tetrahedra[i][3] = d;
            for (unsigned int j=0;j<4;++j)
                if ((unsigned)tetrahedra[i][j] >= points.size())
                {
                    std::cerr << "ERROR: invalid index " << tetrahedra[i][j] << " in tetrahedron " << i << std::endl;
                    tetrahedra[i][j] = 0;
                }
        }
    }

    return true;
}

template<class VecCoord, class VecTetra, class VecTri>
bool read_mesh_vtk(const char* filename, VecCoord& points, VecTetra& tetrahedra, VecTri& triangles)
{
    std::ifstream in(filename);
    if (!in)
    {
        std::cerr << "Cannot open file " << filename << std::endl;
        return false;
    }
    std::cout << "Reading file " << filename << std::endl;
    std::string type;
    while (in >> type) {
        if (type == "POINTS") {
            if (! read_vtk_node(in,points)) return false;
        } else if (type == "CELLS") {
            if (! read_vtk_element(in,points,tetrahedra,triangles)) return false;
        }
    }

    in.close();

    return true;
}


#endif

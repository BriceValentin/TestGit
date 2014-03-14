#include "simulation.h"
#include <DataStructures.h>

#include <fstream>
#include<string>
#include<iostream>
#include <cstring>

bool compare(const char * key,std::string & line) {
    int len = strlen(key);

    return ! strncmp(key,line.c_str(),len);
}

int getAttributOffset(std::string & line) {
    const char * c_line = line.c_str();

    int id=0;
    while (c_line[id] > 0 && c_line[id]!='=') id++;
    id++;
    while (c_line[id] > 0 && c_line[id]==' ') id++;
    //id++;

    return id;
}

std::string getString(std::string & line) {
    return std::string(line.c_str() + getAttributOffset(line));
}

double getDouble(std::string & line) {
    std::istringstream is_param(line.c_str() + getAttributOffset(line));
    double res;
    is_param>>res;
    return res;
}

int getInt(std::string & line) {
    std::istringstream is_param(line.c_str() + getAttributOffset(line));
    int res;
    is_param>>res;
    return res;
}

bool getBool(std::string & line) {
    std::istringstream is_param(line.c_str() + getAttributOffset(line));
    bool res;
    is_param>>res;
    return res;
}

TimeIntegration getIntegration(std::string & line) {
    std::string attr = getString(line);
    if (compare("ODE_EulerExplicit",attr)) {
        return ODE_EulerExplicit;
    } else if (compare("ODE_EulerImplicit",attr)) {
        return ODE_EulerImplicit;
    } else {
        std::cerr << "Unknown integration in " << line << " use ODE_EulerExplicit default ODE_EulerExplicit" << std::endl;
        return ODE_EulerExplicit;
    }
}

TVec3 getVec3(std::string & line) {
    std::istringstream is_param(line.c_str() + getAttributOffset(line));
    double v1;
    double v2;
    double v3;

    is_param>>v1;
    is_param>>v2;
    is_param>>v3;

    TVec3 res(v1,v2,v3);

    return res;
}

BoundingBox getBBox(std::string & line) {
    std::istringstream is_param(line.c_str() + getAttributOffset(line));
    double v1;
    double v2;
    double v3;
    double v4;
    double v5;
    double v6;
    is_param>>v1;
    is_param>>v2;
    is_param>>v3;
    is_param>>v4;
    is_param>>v5;
    is_param>>v6;

    BoundingBox res;
    res.min = TVec3(v1,v2,v3);
    res.max = TVec3(v4,v5,v6);

    return res;
}

bool  read_config_file(std::string & parentDir,std::string & filename)
{
    std::cout << "Read config file " << filename << std::endl;
    std::ifstream file;
    file.open((parentDir + filename).c_str());

    if (! file.is_open()) {
        std::cerr << "Error cannot find the file " << filename << std::endl;
        return false;
    }

    std::string line;
    while (getline(file, line)) {
        if (compare("//",line) || line.length() == 0) {
            //do nothing its a comment
        } else if (compare("msh_filename",line)) {
            simulation_params.fem_msh_filename = parentDir + getString(line);
            simulation_params.fem_vtk_filename.clear();
            simulation_params.fem_netget_filename.clear();
        } else if (compare("fem_filename",line)) {
            simulation_params.fem_netget_filename = parentDir + getString(line);
            simulation_params.fem_vtk_filename.clear();
            simulation_params.fem_msh_filename.clear();
        } else if (compare("vtk_filename",line)) {
            simulation_params.fem_vtk_filename = parentDir + getString(line);
            simulation_params.fem_netget_filename.clear();
            simulation_params.fem_msh_filename.clear();
        } else if (compare("render_filenames",line)) {
            simulation_params.render_filenames.push_back(parentDir + getString(line));
        } else if (compare("colorMap",line)) {
            simulation_params.colorMap = parentDir + getString(line);
        } else if (compare("normalMap",line)) {
            simulation_params.normalMap = parentDir + getString(line);
        } else if (compare("shaderVertex",line)) {
            simulation_params.shaderVertex = parentDir + getString(line);
        } else if (compare("shaderFragment",line)) {
            simulation_params.shaderFragment = parentDir + getString(line);
        } else if (compare("timeStep",line)) {
            simulation_params.timeStep = getDouble(line);
        } else if (compare("subSteps",line)) {
            simulation_params.subSteps = getInt(line);
        } else if (compare("odeSolver",line)) {
            simulation_params.odeSolver = getIntegration(line);
        } else if (compare("rayleighMass",line)) {
            simulation_params.rayleighMass = getDouble(line);
        } else if (compare("rayleighStiffness",line)) {
            simulation_params.rayleighStiffness = getDouble(line);
        } else if (compare("maxIter",line)) {
            simulation_params.maxIter = getInt(line);
        } else if (compare("tolerance",line)) {
            simulation_params.tolerance = getDouble(line);
        } else if (compare("massDensity",line)) {
            simulation_params.massDensity = getDouble(line);
        } else if (compare("youngModulus",line)) {
            simulation_params.youngModulus = getDouble(line);
        } else if (compare("poissonRatio",line)) {
            simulation_params.poissonRatio = getDouble(line);
        } else if (compare("reorder",line)) {
            simulation_params.reorder = getBool(line);
        } else if (compare("flip_tetra",line)) {
            simulation_params.flip_tetra = getBool(line);
        } else if (compare("gravity",line)) {
            simulation_params.gravity = getVec3(line);
        } else if (compare("fixedBBox",line)) {
            simulation_params.fixedBBox.push_back(getBBox(line));
        } else {
            std::cerr << "Unknown parameter " << line << std::endl;
        }
    }

    return true;
}


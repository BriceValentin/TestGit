#include <simulation.h>
#include <rendering/render.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <rendering/glut_methods.h>
#include <read_config.h>

std::string parentDir;

/// Get the full path of the current process. The given filename should be the value of argv[0].
std::string getProcessFullPath(const char* filename)
{
    if (!filename || filename[0]!='/')
    {
        char path[1024];
        memset(path,0,sizeof(path));
        if (readlink("/proc/self/exe",path,sizeof(path)-1) == -1)
          std::cerr <<"Error: can't read the contents of the link." << std::endl;
        if (path[0])
            return path;
        else
            std::cout << "ERROR: can't get current process path..." << std::endl;
    }

    return filename;
}

std::string getParentDir(std::string path)
{
    std::string::size_type pos = path.find_last_of("/\\");
    if (pos == std::string::npos)
        return ""; // no directory
    else
        return path.substr(0,pos);
}

static void cpuid(unsigned int a, unsigned int b[4])
{
    asm volatile("xchgl %%ebx, %1\n"
                 "cpuid\n"
                 "xchgl %%ebx, %1\n"
                 :"=a"(*b),"=r"(*(b+1)),
                 "=c"(*(b+2)),"=d"(*(b+3)):"0"(a));
}

std::string cpu_name()
{
    unsigned int b[13] = {0};
    cpuid(0x80000000,b);
    unsigned int max = b[0];
    if (max < 0x80000004) return std::string();
    cpuid(0x80000002,b);
    cpuid(0x80000003,b+4);
    cpuid(0x80000004,b+8);
    std::string s;
    b[12] = 0;
    const char* p = (const char*)b;
    char last = '\0';
    while (*p)
    {
        char c = *p; ++p;
        if (c == ' ' && last == ' ') continue;
        if (c == '(')
        {
            while (*p && c != ')') c = *p++;
            continue;
        }
        s += c; last = c;
    }
    return s;
}

bool simulation_preload(int argc, char **argv) {
    parentDir = getProcessFullPath(argv[0]);
    parentDir = getParentDir(parentDir);
    if (!parentDir.empty()) parentDir += '/';

    std::ostringstream o;
    o << "CPU: " << cpu_name();
    simulation_params.device_name = o.str();

    simulation_params.timeStep = 0.05;
    simulation_params.subSteps = 1;
    simulation_params.odeSolver = ODE_EulerImplicit;
    simulation_params.rayleighMass = 0.1;
    simulation_params.rayleighStiffness = 0.1;
    simulation_params.maxIter = 50;
    simulation_params.tolerance = 1e-3;
    simulation_params.massDensity = 0.1;
    simulation_params.gravity = TVec3(0,-10,0);
    simulation_params.youngModulus = 500000;
    simulation_params.poissonRatio = 0.4;
    simulation_params.reorder = true;
    simulation_params.flip_tetra = false;

    if (argc>1) {
        std::string config(argv[1]);
        return read_config_file(parentDir,config);
    } else {
        /// Default parameters
        simulation_params.fem_msh_filename = parentDir + "data/liver.msh";
        simulation_params.render_filenames.push_back(parentDir + "data/liver-smoothUV.obj");

        simulation_params.colorMap = parentDir + "data/liver-cyrrhosis.png";
        simulation_params.normalMap = parentDir + "data/liver-cyrrhosis-NM.png";

        simulation_params.shaderVertex = parentDir + "data/bumpVertexShader.glsl";
        simulation_params.shaderFragment = parentDir + "data/bumpFragmentShader.glsl";

        BoundingBox bb;
        bb.min = TVec3(-3,5,-3);
        bb.max = TVec3(3,6,3);
        simulation_params.fixedBBox.push_back(bb);

        return true;
    }
}

bool main_load()
{
    std::cout << "Load meshes" << std::endl;

    if (!simulation_load_fem_mesh()) return false;

    if (!simulation_load_render_mesh()) return false;

    if (simulation_params.reorder) simulation_reorder_fem_mesh();

    std::cout << "Init simulation" << std::endl;

    if (!simulation_init()) return false;

    return true;
}

int main(int argc, char **argv)
{
    if (!simulation_preload(argc,argv)) return 1;

    init_glut(&argc, argv);

    if (! main_load()) return 1;

    run_glut();

    return 0;
}




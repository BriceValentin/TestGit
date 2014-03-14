
TEMPLATE = app
TARGET = FEM
OBJECTS_DIR = ../OBJ
DESTDIR = ../bin

HEADERS = src/mesh/read_mesh_netgen.h \
          src/mesh/read_mesh_obj.h \
          src/mesh/write_mesh_obj.h \
          src/mesh/FEMmesh.h \
          src/mesh/Surfacemesh.h \
          src/rendering/glut_methods.h \
          src/rendering/render.h \
          src/octree/octree.h \
          src/mapping/Mapping_kernels.h \
          src/read_config.h \
          src/simulation.h
	  
SOURCES = src/mesh/FEMmesh.cpp \
          src/mesh/Surfacemesh.cpp \
          src/rendering/glut_methods.cpp \
          src/rendering/gui.cpp \
          src/rendering/render.cpp \
          src/mapping/BarycentricMapping.cpp \
          src/main.cpp \
          src/read_config.cpp \
          src/simulation.cpp

INCLUDEPATH += src
INCLUDEPATH += ../DataStructures/src
INCLUDEPATH += ../State/src
INCLUDEPATH += ../Mass/src
INCLUDEPATH += ../ForceField/src
INCLUDEPATH += ../ExplicitIntegration/src
INCLUDEPATH += ../ImplicitIntegration/src

LIBS += -L../bin/

LIBS += -lglut \
	-lGLU \
	-ldl \
        -lGLEW \
        -lpng \
        -lDataStructures \
        -lState \
        -lMass \
        -lForceField \
        -lExplicitIntegration \
        -lImplicitIntegration

DEPENDPATH = $$INCLUDEPATH



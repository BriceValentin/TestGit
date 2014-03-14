
TEMPLATE = lib
TARGET = ExplicitIntegration
OBJECTS_DIR = ../OBJ
DESTDIR = ../bin

HEADERS = src/ExplicitIntegration_kernels.h
	  
SOURCES = src/ExplicitIntegration.cpp

LIBS += -L../bin/

INCLUDEPATH += src
INCLUDEPATH += ../DataStructures/src
INCLUDEPATH += ../ForceField/src
INCLUDEPATH += ../Mass/src
INCLUDEPATH += ../State/src

LIBS += -lDataStructures \
        -lForceField \
        -lMass \
        -lState \

DEPENDPATH = $$INCLUDEPATH


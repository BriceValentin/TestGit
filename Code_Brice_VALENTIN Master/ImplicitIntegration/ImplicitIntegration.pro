
TEMPLATE = lib
TARGET = ImplicitIntegration
OBJECTS_DIR = ../OBJ
DESTDIR = ../bin

HEADERS = src/ImplicitIntegration_kernels.h
	  
SOURCES = src/ImplicitIntegration.cpp

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

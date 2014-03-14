
TEMPLATE = lib
TARGET = ForceField
OBJECTS_DIR = ../OBJ
DESTDIR = ../bin

HEADERS = src/ForceField_kernels.h
	  
SOURCES = src/TetrahedronFEMForceField.cpp

INCLUDEPATH += src
INCLUDEPATH += ../DataStructures/src

LIBS += -L../bin/

LIBS += -lDataStructures

DEPENDPATH = $$INCLUDEPATH

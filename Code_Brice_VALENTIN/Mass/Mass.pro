
TEMPLATE = lib
TARGET = Mass
OBJECTS_DIR = ../OBJ
DESTDIR = ../bin

HEADERS = src/Mass_kernels.h
	  
SOURCES = src/UniformMass.cpp

INCLUDEPATH += src
INCLUDEPATH += ../DataStructures/src

LIBS += -L../bin/

LIBS += -lDataStructures

DEPENDPATH = $$INCLUDEPATH

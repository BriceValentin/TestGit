
TEMPLATE = lib
TARGET = State
OBJECTS_DIR = ../OBJ
DESTDIR = ../bin

HEADERS = src/State_kernels.h
	  
SOURCES = src/MechanicalObject.cpp \
          src/FixedConstraint.cpp \

INCLUDEPATH += src
INCLUDEPATH += ../DataStructures/src

LIBS += -L../bin/

LIBS += -lDataStructures \
        -lForceField \
        -lMass \

DEPENDPATH = $$INCLUDEPATH

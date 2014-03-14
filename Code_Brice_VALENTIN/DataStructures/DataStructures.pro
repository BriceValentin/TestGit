
TEMPLATE = lib
TARGET = DataStructures
OBJECTS_DIR = ../OBJ
DESTDIR = ../bin

HEADERS = src/DataStructures.h \
          src/Sofa/DataTypeInfo.h \
          src/Sofa/fixed_array.h \
          src/Sofa/Mat.h \
          src/Sofa/Vec.h \
	  
SOURCES = src/DataStructures.cpp

INCLUDEPATH += src

DEPENDPATH = $$INCLUDEPATH

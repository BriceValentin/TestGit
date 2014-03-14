TEMPLATE = subdirs

INCLUDEPATH += Simulation/src
INCLUDEPATH += DataStructures/src
INCLUDEPATH += ForceField/src
INCLUDEPATH += Mass/src
INCLUDEPATH += ExplicitIntegration/src
INCLUDEPATH += ImplicitIntegration/src

SUBDIRS += DataStructures
SUBDIRS += Mass
SUBDIRS += ForceField
SUBDIRS += State
SUBDIRS += ExplicitIntegration
SUBDIRS += ImplicitIntegration
SUBDIRS += Simulation

OTHER_FILES = bin/config

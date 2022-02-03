AUTOMAKE_OPTIONS = foreign
bin_PROGRAMS = crackMesher

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = bootstrap

AM_CPPFLAGS = $(BOOST_CPPFLAGS) -std=c++14

AM_LDFLAGS = $(BOOST_LDFLAGS)

crackMesher_LDADD = \
    $(BOOST_FILESYSTEM_LIB) \
    $(BOOST_SYSTEM_LIB) \
    $(BOOST_PROGRAM_OPTIONS_LIB)
    
crackMesher_SOURCES = \
src/crackMesher.cpp \
src/GlobalMesh.cpp \
src/Grain.cpp \
src/refine/coarsenSurf.h \
src/refine/refineSurf.h \
src/MeshIO/asciiSTL.cpp \
src/MeshIO/INPReader.cpp \
src/MeshIO/TetGen.cpp \
src/programSettings.h

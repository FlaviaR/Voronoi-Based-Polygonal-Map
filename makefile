# For compiling for windows:
#    make OS="win32"
# For cross-compiling for windows:
#    make OS="cross"

ARCH = $(shell getconf LONG_BIT)
OS = linux
CC = c++
WARCH=w64

.SUFFIXES: .C .cpp .cxx .h

# include $(CGAL_MAKEFILE)

LDFLAGS = -lCGAL
CXX     = g++

ifeq ($(OS),linux)
   # flags for C++ compiler:
   CFLAGS  = -I/usr/include/ -frounding-math
   # libraries to link with:
   LIBPATH = $(CGAL_LIBPATH) -L/usr/lib$(ARCH)
   OGLLIBS = -lglut -lGLU -lGL
   EXT = "-linux$(ARCH)"
else ifeq ($(OS),osx)
   # flags for C++ compiler:
   CFLAGS  = -I/opt/local/include/ -frounding-math
   # libraries to link with:
   LIBPATH = $(CGAL_LIBPATH) -L/opt/local/lib -L/usr/X11/lib
   OGLLIBS = -lglut -lGLU -lGL
   EXT=".osx"
else ifeq ($(OS),win32)
   CXX = i686-w64-mingw32-g++
   # flags for C++ compiler:
   CFLAGS  = -frounding-math
   # libraries to link with:
   LIBPATH =
   OGLLIBS = -lglut -lglu32 -lopengl32 -lboost_system
   EXT=".exe"
   OPTIONS=-static-libgcc -static-libstdc++
endif

# **********************************************************************************
all:	VoronoiDemo$(EXT)

VoronoiDemo$(EXT):	VoronoiDemo.o
	$(CXX)  -o $@ $^ $(CFLAGS) $(LIBPATH) $(LDFLAGS) $(OGLLIBS) $(OPTIONS)

.c.o:	$*.h
	$(CXX) -c $(CFLAGS) $*.c

.cpp.o:	$*.h
	$(CXX) -c $(CFLAGS) $*.cpp

clean:
	rm *.o VoronoiDemo$(EXT)

.C.o:	$*.h
	$(CXX) -c $(CFLAGS) $*.C



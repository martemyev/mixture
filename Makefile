CXX=g++ # c++ compiler
CXXFLAGS=-c -Wall -O3 -I/home/artemiev/projects/multFEMv1.0# -w1 -wd981 -DDEBUG
LDFLAGS=
LIBSOURCES=Geometry.cpp TriangularMesh.cpp pugixml.cpp GeoShape.cpp Cylinder.cpp Ellipsoid.cpp OrthoBrick.cpp GeoTetrahedron.cpp
LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBRARY=libmixture.a

all: $(LIBOBJECTS)
	ar -cvq $(LIBRARY) $(LIBOBJECTS)
	ranlib -t

# build object files
.cpp.o:
	$(CXX) $(CXXFLAGS) $^

clean:
	rm -rf $(LIBOBJECTS) $(LIBRARY)

rebuild: clean all

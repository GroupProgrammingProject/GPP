CXX = g++
CXXFLAGS=
LD = g++

#files to be compiled in libraries
LIBFILES = kpointsfunctions MolDyn band_hamiltonian hamiltonian geometryinfo readinxyz vectorfunctions ScaleGeom Gethijab functions 
LIBOBJECTS = $(addsuffix .o, $(LIBFILES))
LDLIBS =  $(addprefix -l,$(LIBFILES))
LIBNAMES =  $(addsuffix .so, $(LIBFILES))

#include path to eigen
INCLUDE = $(addprefix -I , ~/Software/eigen/)
HEADERS = ../include/*.h

#compile and runtime shared lib linking path
LDFLAGS = $(addprefix -L ,../lib/)
RPATH = -Wl,-rpath=../lib 

# runs when make is executed without further options
# meant to compile executable with dynamic lib made by make install
all: genkgrid genkpath replicate scalecell

genkgrid: gen_kgrid.cpp
	$(CXX) $(INCLUDE) $(RPATH) $(LDFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
genkgrid: $(HEADERS)

genkpath: gen_kpath.cpp
	$(CXX) $(INCLUDE) $(RPATH) $(LDFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
genkpath: $(HEADERS)

replicate: ReplicateCell.cpp 
	$(CXX) $(INCLUDE) $(RPATH) $(LDFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
replicate: $(HEADERS)

scalecell: ScaleCell.cpp
	$(CXX) $(INCLUDE) $(RPATH) $(LDFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
scalecell: $(HEADERS)

# install, needs to run to compile shared library objects
#install: $(LIBNAMES) 

#%.so: ./%.o
#	$(CXX) -shared -fPIC -o lib/lib$@ $<

%.o: ../src/%.cpp
	$(CXX) -c $(INCLUDE) -fPIC $<
%.o: ../include/%.h

.PRECIOUS: %.o 

clean: 
	rm *.o
CXX = g++
CXXFLAGS= 
LD = g++

#files to be compiled in libraries
LIBFILES = MolDyn hamiltonian geometryinfo readinxyz vectorfunctions ScaleGeom Gethijab functions 
LIBOBJECTS = $(addsuffix .o, $(LIBFILES))
LDLIBS =  $(addprefix -l,$(LIBFILES))
LIBNAMES =  $(addsuffix .so, $(LIBFILES))
LIBPATH = ./lib

#include path to eigen
HEADERPATH = ./include/
INCLUDE = $(addprefix -I , ~/Software/eigen/ $(HEADERPATH))
HEADERS = ./include/*.h

#compile and runtime shared lib linking path
LDFLAGS = $(addprefix -L ,./lib/)
RPATH = -Wl,-rpath=./lib 

# runs when make is executed without further options
# meant to compile executable with dynamic lib made by make install
all: main

main: relax_main.cpp
	$(CXX) $(LDFLAGS) $(INCLUDE) $(RPATH) $(CXXFLAGS) $< -o $@ $(LDLIBS)
main: $(HEADERS) $(LIBPATH)/*.so


# install, needs to run to compile shared library objects
install: $(LIBNAMES) 

%.so: %.o
	$(CXX) -shared -fPIC -o $(LIBPATH)/lib$@ $<

%.o: src/%.cpp
	$(CXX) -c $(INCLUDE) -fPIC $<
#%.o: include/%.h

clean: 
	rm *.o

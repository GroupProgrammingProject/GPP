CXX = g++
CXXFLAGS=
LD = g++

#files to be compiled in libraries
LIBFILES = functions geometryinfo hamiltonian MolDyn readinxyz vectorfunctions ScaleGeom Gethijab
LIBOBJECTS = $(addsuffix .o, $(LIBFILES))
LDLIBS =  $(addprefix -l,$(LIBFILES))
LIBNAMES =  $(addsuffix .so, $(LIBFILES))

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

main: main.cpp
	$(CXX) $(INCLUDE) $(RPATH) $(LDFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)
main: $(HEADERS)


# install, needs to run to compile shared library objects
install: $(LIBNAMES) 

%.so: %.o
	$(CXX) -shared -fPIC -o lib/lib$@ $<

%.o: src/%.cpp
	$(CXX) -c $(INCLUDE) -fPIC $<
%.o: include/%.h

.PRECIOUS: %.o 

clean: 
	rm *.o

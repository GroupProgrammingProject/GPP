CXX = g++
CXXFLAGS+= -std=c++11
LD = g++

#files to be compiled in libraries
LIBFILES = phonons kpointsfunctions MolDyn band_hamiltonian hamiltonian geometryinfo readinxyz readinxyzv vectorfunctions ScaleGeom Gethijab functions 
LIBOBJECTS = $(addsuffix .o, $(LIBFILES))
LDLIBS =  $(addprefix -l,$(LIBFILES))
LIBNAMES =  $(addsuffix .so, $(LIBFILES))

# Directories that will be created
LIB_DIR = ./lib

#include path to eigen
INCLUDE = $(addprefix -I , ~/Software/eigen/)
HEADERS = ./include/*.h

#compile and runtime shared lib linking path
LDFLAGS = $(addprefix -L ,$(LIB_DIR))
RPATH = -Wl,-rpath=$(LIB_DIR) 

# runs when make is executed without further options
# meant to compile executable with dynamic lib made by make install
all: singleE_main md_main phonons_main relax_main

singleE_main: singleE_main.cpp
	$(CXX) $(LDFLAGS) $(INCLUDE) $(RPATH) $(CXXFLAGS) $< -o $@ $(LDLIBS)
singleE_main: $(HEADERS)

md_main: md_main.cpp
	$(CXX) $(LDFLAGS) $(INCLUDE) $(RPATH) $(CXXFLAGS) $< -o $@ $(LDLIBS)
md_main: $(HEADERS)

phonons_main: phonons_main.cpp
	$(CXX) $(LDFLAGS) $(INCLUDE) $(RPATH) $(CXXFLAGS) $< -o $@ $(LDLIBS)
phonons_main: $(HEADERS)

relax_main: relax_main.cpp
	$(CXX) $(LDFLAGS) $(INCLUDE) $(RPATH) $(CXXFLAGS) $< -o $@ $(LDLIBS)
relax_main: $(HEADERS)

# install, needs to run to compile shared library objects
install: $(LIB_DIR) $(LIBNAMES) 

$(LIB_DIR): 
	mkdir -p $@

%.so: %.o
	$(CXX) -shared -fPIC -o lib/lib$@ $<

%.o: src/%.cpp
	$(CXX) -c $(INCLUDE) -fPIC $< 
%.o: include/%.h

.PRECIOUS: %.o 

clean: 
	rm *.o *.txt

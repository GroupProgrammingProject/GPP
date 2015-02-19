CXX = g++
CXXFLAGS=$(CFLAGS)
LD = g++
INCLUDE = $(addprefix -I ,~/Software/eigen/)
HEADERS = ./include/*.h

main.o:main.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $(CXXFLAGS) $< -o $@
main.o: $(HEADERS)

main: main.o
	$(LD) -o $@ $(LDFLAGS) -I  main.o

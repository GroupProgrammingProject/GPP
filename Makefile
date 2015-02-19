CXX = g++
CXXFLAGS=$(CFLAGS)
LD = g++
INC = -I ~/Software/eigen/
HEADERS = ./include/*.h

main: main.cpp
	$(LD) $(LDFLAGS) $(INC) -o $@ $^
main.o: $(HEADERS)

CXX = g++
LD = g++
scaling_functions: functions.o main.o
	 $(LD) $(LDFLAGS) -o $@ $^

main.o:functions.h

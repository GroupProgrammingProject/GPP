#ifndef VECTORFUNCTIONS_H
#define VECTORFUNCTIONS_H

// Some vector functions not included in the std::vector library
#include <iostream>
#include <fstream>
#include <vector>

template <typename T> void Print(std::vector<T>* vect);
template <typename T> void Print(std::vector<T>* vect, char* ofilename);

#endif

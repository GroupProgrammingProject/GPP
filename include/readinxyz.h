#ifndef READINXYZ_H
#define READINXYZ_H

// Read in from xyz file 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

void ReadInXYZ (char* filename, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* lats, bool pbc);
std::map<std::string, int> SetElementMap();

#endif

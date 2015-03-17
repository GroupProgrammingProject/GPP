#ifndef READINXYZV_H
#define READINXYZV_H

// Read in from xyz file 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <stdlib.h>
#include <sstream>
void pbcshift(std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* lats);

void ReadInXYZV (char* filename, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz,std::vector<double>* lats, bool pbc, std::vector<bool>* velspec);
std::map<std::string, int> SetElementMap();

#endif

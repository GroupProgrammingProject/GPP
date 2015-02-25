#ifndef READINXYZ_H
#define READINXYZ_H

// Read in from xyz file 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

// Reads in 4 vectors from cell file: elements, x , y and z coords
void ReadInXYZ (char* filename, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z);
// Function to initiate a map to translate element symbol to integer values corresponding to its atomic number
std::map<std::string, int> SetElementMap();

#endif

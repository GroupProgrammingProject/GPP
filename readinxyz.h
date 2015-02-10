#ifndef READINXYZ_H
#define READINXYZ_H

// Read in from xyz file 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

void ReadInXYZ (char* filename, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z);
std::map<std::string, int> SetElementMap();

// Function to initiate a map to translate element symbol to integer values corresponding to its atomic number
std::map<std::string, int> SetElementMap(){
	std::map<std::string, int> map;
	map["H"]=1;
	map["C"]=6;
	return map;
}

// Reads in 4 vectors from cell file: elements, x , y and z coords
void ReadInXYZ(char* filename,std::vector<double>* xvect, std::vector<double>* yvect, std::vector<double>* zvect){
	std::ifstream infile(filename);
	// Initiate the map of elements (symbol-> atomic number)
	std::map<std::string, int> elmap=SetElementMap();
	// First two lines are handled seperately and are not to be stored
	int num;
	char skip[2];
	infile>>num;
	infile>>skip>>skip>>skip>>skip;
	// Store the molecule type and x, y, z positions
	double x, y, z;
	std::string type;
	while (infile>>type>>x>>y>>z){
		xvect->push_back(x);
		yvect->push_back(y);
		zvect->push_back(z);
	}
}
#endif

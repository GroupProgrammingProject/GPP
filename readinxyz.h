// Read in from file functionality
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Reads in 4 vectors from cell file
void ReadInXYZ(char* filename, std::vector<std::string>* names, std::vector<double> x, std::vector<double> y, std::vector<double> z){
	std::ifstream infile;
	infile.open(filename);
	std::string textarray[10];
}

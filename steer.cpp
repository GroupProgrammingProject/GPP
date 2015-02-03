// File that assembles all the different TB and MD modules
#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "readinxyz.h"
#include "vectorfunctions.h"
#include "hamiltonian.h"
#include "functions.h"
#include "Gethijab.h"
#include "MolDyn.h"

int main(int argc, char* argv[]){
	// Read in types, 
	std::vector<int> types;
	std::vector<double> X, Y, Z;
	std::vector<int>* typeptr = &types;
	std::vector<double>* Xptr = &X;
	std::vector<double>* Yptr = &Y;
	std::vector<double>* Zptr = &Z;
	ReadInXYZ (argv[1], typeptr, Xptr, Yptr, Zptr);
// Debug output
Print (typeptr);
Print (Xptr);
Print (Yptr);
Print (Zptr);
	
return 0;
}

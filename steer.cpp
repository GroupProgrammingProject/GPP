// File that assembles all the different TB and MD modules
#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "readinxyz.h"
#include "vectorfunctions.h"
#include "hamiltonian.h"
#include "functions.h"
//#include "MolDyn.h"

int main(int argc, char* argv[]){
	// Read in types, 
	std::vector<int> type;
	std::vector<double> posx, posy, posz;
	std::vector<int>* typeptr = &type;
	std::vector<double>* Xptr = &posx;
	std::vector<double>* Yptr = &posy;
	std::vector<double>* Zptr = &posz;
	ReadInXYZ (argv[1], typeptr, Xptr, Yptr, Zptr);
// Debug output
/*Print (typeptr);
Print (Xptr);
Print (Yptr);
Print (Zptr);
*/
	double ebs;
	// Number of atoms
	int n=type.size();
	Matrix Hamiltonian;
	Matrix* Hamiltonianptr=&Hamiltonian;

	ebs=Hamiltonian(n,typeptr,Xptr,Yptr,Zptr);
	std::cout << "Ebs=" << ebs << " eV" << std::endl;
	
return 0;
}

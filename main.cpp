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
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	// Read in types, 
	std::vector<int> type;
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz);
	// Number of atoms
	int n=type.size();
	// Create empty eigenvector array needed for MD
	std::vector<double> eigvects(16*n*n);
	// Energies from TB model
	double ebs,erep,etot;

	// Starting TB	module: calculating energies and eigenvectors
	ebs=Hamiltonian(n,&type,&posx,&posy,&posz,&eigvects);
	erep=Erep(&type,&posx,&posy,&posz);
// Determining erep works, however, we need to check if we're passing pointers or arrays to Erep()
	etot=ebs+erep;

return 0;
}

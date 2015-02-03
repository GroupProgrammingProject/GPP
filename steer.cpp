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
	ReadInXYZ (argv[1], &type, &posx, &posy, &posz);
	// Number of atoms
	int n=type.size();
	// Create empty arrays needed for MD
	std::vector<double> H_MD(16*n*n),eigvects(16*n*n);
	// Energies from TB model
	double ebs,erep,etot;

	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&type,&posx,&posy,&posz,&H_MD,&eigvects);
	//H_MD and eigvects have now also been populated
	erep=Erep(type,posx,posy,posz);
// Determining erep works, however, we need to check if we're passing pointers or arrays to Erep()
	etot=ebs+erep;
	

return 0;
}

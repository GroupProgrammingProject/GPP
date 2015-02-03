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
	std::vector<double> eigvects(16*n*n);
	// Energies from TB model
	double ebs,erep,etot;

	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&type,&posx,&posy,&posz,&eigvects);
	//H_MD and eigvects have now also been populated
	erep=Erep(type,posx,posy,posz);
// Determining erep works, however, we need to check if we're passing pointers or arrays to Erep()
	etot=ebs+erep;

//test energy
	std::cout << "Ebs = " << ebs << std::endl;

	// Start MD routine: steepest descent
	
	// First call near_neigh()
	int nnmax=10;	
	std::vector<double> mass(n),nnear(n),inear(n*nnmax);		// Mass in kg?
	double rc=2.6;
	double sx=10,sy=10,sz=10;
//	near_neigh(n,&posx,&posy,&posz,rc,&nnear,&inear,sx,sy,sz);

	// Start forces convergence loop
	int norbs=4;
	std::vector<double> fx(n),fy(n),fz(n);
/*	
	while(convergence limit)
	forces(n,norbs,&posx,&posy,&posz,&eigvects,rc,&nnear,&inear,&fx,&fy,&fz);
	change positions of atoms

*/

	//Plot positions!
	

return 0;
}

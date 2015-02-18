#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "readinxyz.h"
#include "vectorfunctions.h"
#include "hamiltonian.h"
#include "functions.h"
#include "geometryinfo.h"
#include "ScaleGeom.h"

int main(int argc, char* argv[]){
	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0;
	// Read in types, 
	std::vector<int> type;
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &type, &posx, &posy, &posz);
	// Number of atoms
	int n=posx.size();
	std::cout << "n = " << n <<std::endl;
	std::vector<double> posxnew(n), posynew(n), posznew(n);

	// Scale cell size
	//ScaleGeom(n, 1.2, 1.2, 1.2, &posx, &posy, &posz, &posxnew, &posynew, &posznew);

	// Calculate distances
	std::vector<double> modr(n*n),rx(n*n),ry(n*n),rz(n*n);
	GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz);
	// Create empty arrays needed for MD
	std::vector<double> eigvects(16*n*n);
	// Energies from TB model
	double ebs,erep,etot;

	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
// Determining erep works, however, we need to check if we're passing pointers or arrays to Erep()
	etot=ebs+erep;

	std::cout << "Ebs = " << ebs << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	std::cout << "Etot = " << etot << std::endl;

return 0;
}

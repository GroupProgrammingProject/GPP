#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "include/readinxyz.h"
#include "include/vectorfunctions.h"
#include "include/hamiltonian.h"
#include "include/functions.h"
#include "include/geometryinfo.h"

int main(int argc, char* argv[]){
	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0;
	// Read in types, 
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz);
	// Number of atoms
	int n=posx.size();
	// Calculate distances
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);

	double a=7.37854,b=8.52,c=2*3.45,rc=2.6,rv=3;

std::cout << "n=" << n << std::endl;

//	GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz);
	PbcGetAllDistances(&modr, &rx, &ry, &rz, &posx, &posy, &posz, a, b, c, rv);

std::cout << "rx" << std::endl;
std::cout << rx << std::endl << std::endl;
std::cout << "ry" << std::endl;
std::cout << ry << std::endl << std::endl;
std::cout << "rz" << std::endl;
std::cout << rz << std::endl << std::endl;
std::cout << "modr" << std::endl;
std::cout << modr << std::endl;

	// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot;

/*
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
	etot=ebs+erep;

	std::cout << "Ebs = " << ebs/(double)n << std::endl;
	std::cout << "Erep = " << erep/(double)n << std::endl;
	std::cout << "Etot = " << etot/(double)n << std::endl;
*/

return 0;
}

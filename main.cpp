#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "readinxyz.h"
#include "vectorfunctions.h"
#include "hamiltonian.h"
#include "functions.h"
#include "geometryinfo.h"

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
	// Calculate distances
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);

	double a=6,b=6,c=6,rc=2.6,rv=3;

std::cout << "All distances" << std::endl;
	GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz);
std::cout << modr << std::endl << std::endl;
/*
std::cout << "Zerod distances >rc" << std::endl;
	GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,rc);
std::cout << modr << std::endl << std::endl;

std::cout << "PBCs:" << std::endl;
	PbcGetAllDistances(&modr, &rx, &ry, &rz, &posx, &posy, &posz, a, b, c, rv);
std::cout << modr << std::endl;
*/
	// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot;


	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
	etot=ebs+erep;

	std::cout << "Ebs = " << ebs << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	std::cout << "Etot = " << etot << std::endl;


return 0;
}

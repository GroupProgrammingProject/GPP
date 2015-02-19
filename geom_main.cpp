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
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz);
	// Number of atoms
	int n=posx.size();
	std::cout << "n = " <<n <<std::endl;
	std::vector<double> posxnew(n), posynew(n), posznew(n);
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	double ebs,erep,etot;

	// Change if using PBCs
	double a =7.0993;
	double b =7.0993;
	double c =7.0993;
	double rv = 2.98;

	Eigen::MatrixXd eigvects(4*n,4*n);

	// Calculate distances
	PbcGetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,a,b,c,rv);
		  
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

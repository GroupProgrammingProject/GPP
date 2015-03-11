#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "include/readinxyz.h"
#include "include/vectorfunctions.h"
#include "include/hamiltonian.h"
#include "include/functions.h"
#include "include/geometryinfo.h"

int main(int argc, char* argv[]){
	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0, renn;
	// Read in types, 
	std::vector<double> lats(3), posx, posy, posz;
	bool pbc = 1;
	ReadInXYZ (argv[1], &posx, &posy, &posz, &lats, pbc);
	// Number of atoms
	int n=posx.size();
	// Calculate distances
	Eigen::MatrixXd modr(n,n), rx(n,n), ry(n,n), rz(n,n);
	GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz);
	// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot;
	// Tight binding parameters
	std::vector<double> TBparam(6);
	// In final version, TB params will be read in from input file
	TBparam[0]=-2.99;		// E_s
	TBparam[1]=3.71;		// E_p
	TBparam[2]=-5;			// V_ss_sigma
	TBparam[3]=4.7;		//	V_sp_sigma
	TBparam[4]=5.5;		// V_pp_sigma
	TBparam[5]=-1.55;		// V_pp_pi

	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&TBparam,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
	etot=ebs+erep;

	std::cout << "Ebs = " << ebs/(double)n << std::endl;
	std::cout << "Erep = " << erep/(double)n << std::endl;
	std::cout << "Etot = " << etot/(double)n << std::endl;

return 0;
}

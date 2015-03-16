#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "include/readinxyz.h"
#include "include/vectorfunctions.h"
#include "include/band_hamiltonian.h"
#include "include/hamiltonian.h"
#include "include/functions.h"
#include "include/geometryinfo.h"
#include "include/kpointsfunctions.h"
#include "include/MolDyn.h"

// Example main for a total energy calculation which uses kpoints

int main(int argc, char* argv[]){
	if (argc<2){std::cout<<"You should append two files to the main object!"<<std::endl;}
	if (argc!=3){std::cout<<"You should append one xyz and one .kpts file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=1;
	std::vector<double> posx, posy, posz;
	std::vector<double> lats(3);
	bool pbc = 1;
	ReadInXYZ (argv[1], &posx, &posy, &posz, &lats, pbc);
	scramble(&posx, &posy, &posz);
	// Number of atoms
	int n=posx.size();
	//std::cout << "n = " <<n <<std::endl;
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	double ebs = 0,erep = 0,etot = 0;
	double rv = 2.98, rc = 2.18;

	Eigen::MatrixXcd eigvects(4*n,4*n);

	Eigen::MatrixXd eigvectsr(4*n,4*n);

	// Calculate distances
	Eigen::MatrixXi inear(n,n);
	std::vector<int> nnear(n);
	std::vector<double> fx(n);
	std::vector<double> fy(n);
	std::vector<double> fz(n);
	std::vector<double> fmag(n);
	std::vector<double> fxtemp(n);
	std::vector<double> fytemp(n);
	std::vector<double> fztemp(n);
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	NearestNeighbours(&inear, &nnear, &modr,rv);

	// Starting TB	module: calculating energies
	std::vector<double> eigenvalaar (4*n);

	std::vector<std::vector<double> > kpoints; 
	readinkpoints(argv[2],&kpoints);
	int ktot = kpoints.size();
	std::cout << "ktot = " << ktot << std::endl;
	std::vector<double> kvec(3);
	 
	double ebstemp;
	for (int i=0;i<ktot;i++) {
	  kvec = kpoints.at(i);
	  std::cout << "\nk = [" << kvec.at(0) << " " << kvec.at(1) << " " << kvec.at(2) << "] " << std::endl;
	  ebstemp=band_Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,&eigenvalaar,&kvec,v);
	  ebs = ebs + ebstemp;
	  std::cout << "Ebs = " << ebstemp << std::endl;
	}

	erep=Erep(&modr);

	std::cout << "Ebs = " << ebs/(double)ktot << std::endl;
	//std::cout << "Erep = " << erep << std::endl;
	etot=ebs/(double)ktot+erep;
	std::cout << "Etot per atom = " << etot/((double)n) << std::endl;

	return 0;
}


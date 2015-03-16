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
#include "include/readeigenmatrix.h"

// Example main for a total energy calculation which uses kpoints

int main(int argc, char* argv[]){
	if (argc<4){std::cout<<"You should append four files to the main object!"<<std::endl;}
	if (argc!=5){std::cout<<"You should append one xyz, one .kpts file and two eigenvector files to the main!!"<<std::endl;}
	
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

	Eigen::MatrixXd reigvects(4*n,4*n);
	reigvects = readMatrix(argv[3]);
	Eigen::MatrixXd ieigvects(4*n,4*n);
	ieigvects = readMatrix(argv[4]);
	Eigen::MatrixXcd eigvects(4*n,4*n);
	//	eigvects = reigvects + std::complex<double>(0,1)*ieigvects;
	eigvects.real() = reigvects;
	eigvects.imag() = ieigvects;
	std::cout << "cvectors = \n" << eigvects << std::endl;

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
	std::cout << "\nrx : \n" << std::endl;
	std::cout << rx << std::endl;
	std::cout << "\nry : \n" << std::endl;
	std::cout << ry << std::endl;
	std::cout << "\nrz : \n" << std::endl;
	std::cout << rz << std::endl;
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
	  //ebstemp=band_Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,&eigenvalaar,&kvec,v);
	  //ebs = ebs + ebstemp;
	  //std::cout << "Ebs = " << ebstemp << std::endl;
	  kforces(n,4,rc, &rx, &ry, &rz, &modr, &eigvects, &nnear, &inear, &fxtemp, &fytemp, &fztemp, &kvec);
	  // Output forces from each kpoint calculation
	  std::cout << "\n Forces " << std::endl;
	  for (int l=0;l<n;l++) {
		 std::cout << l << " " << fxtemp.at(l) << " " << fytemp.at(l) << " " << fztemp.at(l) << std::endl;
	  }
	  // Add forces from different kpoints
	  for (int l=0;l<n;l++) {
		 fx.at(l) = fx.at(l) + (1/(double)ktot)*fxtemp.at(l);
		 fy.at(l) = fy.at(l) + (1/(double)ktot)*fytemp.at(l);
		 fz.at(l) = fz.at(l) + (1/(double)ktot)*fztemp.at(l);
	  }
	}

	for (int i=0;i<n;i++) {
	  fmag.at(i) = sqrt(pow(fx.at(i),2) + pow(fy.at(i),2) + pow(fz.at(i),2));
	}		

	std::cout << "\n Forces " << std::endl;
	for (int i=0;i<n;i++) {
	  std::cout << i << " " << fx.at(i) << " " << fy.at(i) << " " << fz.at(i) << " " << fmag.at(i) << std::endl;
	}	

	//eigvectsr = eigvects.real();
	//forces(n,4,rc, &rx, &ry, &rz, &modr, &eigvectsr, &nnear, &inear, &fxtemp, &fytemp, &fztemp);

	/*	std::cout << "\n Regular Forces " << std::endl;
	for (int i=0;i<n;i++) {
	  std::cout << i << " " << fxtemp.at(i) << " " << fytemp.at(i) << " " << fztemp.at(i) << " " << fmag.at(i) << std::endl;
	  }*/

	return 0;
}


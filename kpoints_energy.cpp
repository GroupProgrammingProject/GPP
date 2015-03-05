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
	bool v=0;
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
	Eigen::MatrixXd eigvectreal(4*n,4*n);
	Eigen::MatrixXcd gammapt(4*n,4*n);
	Eigen::MatrixXcd eigvectstot = Eigen::MatrixXcd::Zero(4*n,4*n);

	// Calculate distances
	Eigen::MatrixXi inear(n,n);
	std::vector<int> nnear(n);
	std::vector<std::complex<double> > fx(n,0.0);
	std::vector<std::complex<double> > fy(n,0.0);
	std::vector<std::complex<double> > fz(n,0.0);
	std::vector<double> fmag(n);
	std::vector<std::complex<double> > fxtemp(n);
	std::vector<std::complex<double> > fytemp(n);
	std::vector<std::complex<double> > fztemp(n);
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
	  kforces(n,4,rc, &rx, &ry, &rz, &modr, &eigvects, &nnear, &inear, &fxtemp, &fytemp, &fztemp, &kvec);
	  // Output forces from each kpoint calculation
	  std::cout << "\n Forces " << std::endl;
	  for (int l=0;l<n;l++) {
		 std::cout << l << " " << fxtemp.at(l) << " " << fytemp.at(l) << " " << fztemp.at(l) << std::endl;
	  }
	  // Output real components of forces from kpoints
	  /*std::cout << "\n Real forces:" << std::endl;
	  for (int l=0;l<n;l++) {
		 std::cout << l << " " << fxtemp.at(l).real() << " " << fytemp.at(l).real() << " " << fztemp.at(l).real() << std::endl;
		 }*/
	  // Add forces from different kpoints
	  for (int l=0;l<n;l++) {
		 fx.at(l) = fx.at(l) + std::complex<double>(1/(double)ktot,0)*fxtemp.at(l);
		 //std::cout << "fx.at(l) = " << fx.at(l) << " (1/ktot)fxtemp.at(l) = " << std::complex<double>(1/(double)ktot,0)*fxtemp.at(l) <<std::endl;
		 fy.at(l) = fy.at(l) + std::complex<double>(1/(double)ktot,0)*fytemp.at(l);
		 fz.at(l) = fz.at(l) + std::complex<double>(1/(double)ktot,0)*fztemp.at(l);
	  }
	  //	  std::cout << "\neigenvectors = \n" << eigvects.real() << std::endl;
	  eigvectstot = eigvectstot + (1.0/(double)ktot)*eigvects;
	}

	erep=Erep(&modr);

	std::cout << "Ebs = " << ebs/(double)ktot << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	etot=ebs/(double)ktot+erep;
	std::cout << "Etot per atom = " << etot/((double)n) << std::endl;

	//	std::cout << "\ntotal eigenvectors = \n" << eigvectstot.imag() << std::endl;

	eigvectreal = eigvectstot.real();

	//kforces(n,4,rc, &rx, &ry, &rz, &modr, &eigvectreal, &nnear, &inear, &fx, &fy, &fz);
	
	for (int i=0;i<n;i++) {
	  fmag.at(i) = sqrt(pow(fx.at(i).real(),2) + pow(fx.at(i).imag(),2) + pow(fy.at(i).real(),2) + pow(fy.at(i).imag(),2) + pow(fz.at(i).real(),2) + pow(fz.at(i).imag(),2));
	}		

	std::cout << "\n Forces " << std::endl;
	for (int i=0;i<n;i++) {
	  std::cout << i << " " << fx.at(i) << " " << fy.at(i) << " " << fz.at(i) << " " << fmag.at(i) << std::endl;
	}	
	return 0;
}


#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "include/readinxyz.h"
#include "include/vectorfunctions.h"
#include "include/hamiltonian.h"
#include "include/band_hamiltonian.h"
#include "include/functions.h"
#include "include/geometryinfo.h"
#include "include/kpointsfunctions.h"

int main(int argc, char* argv[]){

	//Turn kpoint sampling on or off
	bool kpts = 0;
	if(kpts==0){
		if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
		if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	}
	if(kpts==1){
		if (argc<2){std::cout<<"You should append two files to the main object!"<<std::endl;}
		if (argc!=3){std::cout<<"You should append one xyz and one .kpts file to the main!!"<<std::endl;}
	}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0;
	// Read in geometry of the structure
	bool pbc=1;
	std::vector<double> lats(3), posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz, &lats, pbc);
	int n=posx.size();
	// Calculate distances
	Eigen::MatrixXd modr(n,n), rx(n,n), ry(n,n), rz(n,n);
	double rv = 2.98;
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create data structures for energy calculations 
	Eigen::MatrixXd eigvects(4*n,4*n);
	double ebs=0,erep=0,etot=0;
	std::vector<double> TBparam(6);

	// In final version, TB params will be read in from input file
	TBparam[0]=-2.99;		// E_s
	TBparam[1]=3.71;		// E_p
	TBparam[2]=-5;			// V_ss_sigma
	TBparam[3]=4.7;		//	V_sp_sigma
	TBparam[4]=5.5;		// V_pp_sigma
	TBparam[5]=-1.55;		// V_pp_pi

	if(kpts==0){
		Eigen::MatrixXd eigvects(4*n,4*n);						//real matrix
		ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
		erep=Erep(&modr);
		etot=ebs+erep;
		std::cout << "Ebs = " << ebs/(double)n << std::endl;
		std::cout << "Erep = " << erep/(double)n << std::endl;
		std::cout << "Etot = " << etot/(double)n << std::endl;
	}
	else if(kpts==1){
		Eigen::MatrixXcd eigvects(4*n,4*n);					//complex matrix
		std::vector<double> eigenvalaar (4*n);
		std::vector<std::vector<double> > kpoints; 
		readinkpoints(argv[2],&kpoints);
		int ktot = kpoints.size();
		std::vector<double> kvec(3);
	 
		int j = 0;
		for (int i=0;i<ktot;i++) {
			kvec = kpoints.at(i);
			ebs=ebs+(1.0/(double)ktot)*band_Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,&eigenvalaar,&kvec,v);
			j++;
		}
		erep=Erep(&modr);
		etot=ebs+erep;
		std::cout << "Ebs = " << ebs/(double)n << std::endl;
		std::cout << "Erep = " << erep/(double)n << std::endl;
		std::cout << "Etot = " << etot/(double)n << std::endl;
	}

return 0;
}

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
	// Read in parameter values passed by the run script with MODE="energy"
	// xyzfilepath=argv[1], kptsfilepath=argv[3]
	int norbs=atoi(argv[8]);
	double rv=atof(argv[6]), rc=atof(argv[7]);
	bool kpts=atoi(argv[2]), pbc=atoi(argv[4]), v=atoi(argv[5]);

	std::vector<double> TBparam(6);
	TBparam[0]=atof(argv[9]);		// E_s
	TBparam[1]=atof(argv[10]);		// E_p
	TBparam[2]=atof(argv[11]);			// V_ss_sigma
	TBparam[3]=atof(argv[12]);		//	V_sp_sigma
	TBparam[4]=atof(argv[13]);		// V_pp_sigma
	TBparam[5]=atof(argv[14]);		// V_pp_pi

/*
	//Turn kpoint sampling on or off
	if(kpts==0){
		if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
		if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	}
	if(kpts==1){
		if (argc<2){std::cout<<"You should append two files to the main object!"<<std::endl;}
		if (argc!=3){std::cout<<"You should append one xyz and one .kpts file to the main!!"<<std::endl;}
	}
*/	
	std::vector<double> lats(3), posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz, &lats, pbc);
	int n=posx.size();
	// Calculate distances
	Eigen::MatrixXd modr(n,n), rx(n,n), ry(n,n), rz(n,n);
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create data structures for energy calculations 
	Eigen::MatrixXd eigvects(4*n,4*n);
	double ebs=0,erep=0,etot=0;
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

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
   //RUNCOMMAND="$0MAIN $1XYZ_FILE_PATH $2KPTS $3KPTS_FILE_PATH $4KSYMM $5PBC $6VERBOSE $7RV $8RC $9NUM_ORBS $10MASS $11P1 $12P2 $13P3 $14P4 $15P5 $16P6"
   int norbs=atoi(argv[9]), maxnn=atoi(argv[10]);
	double rv=atof(argv[7]), rc=atof(argv[8]);
	bool kpts=atoi(argv[2]), pbc=atoi(argv[5]), v=atoi(argv[6]), ksymm=atoi(argv[4]);
	std::vector<double> TBparam(6);
	TBparam[0]=atof(argv[11]);		// E_s
	TBparam[1]=atof(argv[12]);		// E_p
	TBparam[2]=atof(argv[13]);		// V_ss_sigma
	TBparam[3]=atof(argv[14]);		//	V_sp_sigma
	TBparam[4]=atof(argv[15]);		// V_pp_sigma
	TBparam[5]=atof(argv[16]);		// V_pp_pi
	std::vector<double> lats(3), posx, posy, posz, vxin, vyin, vzin;
	std::vector<bool> velspec;
	ReadInXYZ (argv[1], &posx, &posy, &posz, &vxin, &vyin, &vzin, &lats, pbc, &velspec);
	int n=posx.size();
	// Determine maximum number of nearest neighbours
	if(n<maxnn){maxnn=n;}
	// Matrix neighbour list
	Eigen::MatrixXi inear(n,maxnn);
	// Calculate distances
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create data structures for energy calculations 
	Eigen::MatrixXd eigvects(n*norbs,n*norbs);
	double ebs=0,erep=0,etot=0;
	if(kpts==0){
		Eigen::MatrixXd eigvects(n*norbs,n*norbs);						//real matrix
		ebs=Hamiltonian(n,norbs,&TBparam,&modr,&rx,&ry,&rz,&eigvects,v);
	}
	else if(kpts==1){
	  Eigen::MatrixXcd eigvects(n*norbs,n*norbs);					//complex matrix
	  std::vector<double> eigenvalaar (n*norbs);
	  std::vector<std::pair<std::vector<double>,double> > kpoints; 
	  readinkpoints(argv[3],&kpoints,ksymm);
	  std::cout << "ktot = " << kpoints.size() << std::endl;
	  for (int i=0;i<kpoints.size();i++) {
		 std::cout << kpoints.at(i).first.at(0) << "\t" << kpoints.at(i).first.at(1) << "\t" << kpoints.at(i).first.at(2) << "\t" << kpoints.at(i).second << std::endl;
	  }
	  ebs=avekenergy(n,norbs,&rx,&ry,&rz,&modr,&kpoints,&TBparam);
	}
	erep=Erep(&modr);
	etot=ebs+erep;
	std::cout << "Ebs = " << ebs/(double)n << std::endl;
	std::cout << "Erep = " << erep/(double)n << std::endl;
	std::cout << "Etot = " << etot/(double)n << std::endl;
	return 0;
}

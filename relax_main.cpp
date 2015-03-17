#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "include/readinxyz.h"
#include "include/vectorfunctions.h"
#include "include/hamiltonian.h"
#include "include/functions.h"
#include "include/geometryinfo.h"
#include "include/MolDyn.h"

int main(int argc, char* argv[]){
  	// Read in parameter values passed by the run script with MODE="energy"
	// xyzfilepath=argv[1], kptsfilepath=argv[3]
  //RUNCOMMAND="$0MAIN $1XYZ_FILE_PATH $2KPTS $3KPTS_FILE_PATH $4NUM_STEPS $5DT $6PBC $7T $8FRAME_RATE $9VERBOSE $10RV $11RC $12NUM_ORBS $13MAX_NEIGHBOURS $14MASS $15TOL $16MAX_STEEP $17NU $18H $19P1 $20P2 $21P3 $22P4 $23P5 $24P6"
  
	// From now on is the good version (input)
  int nmd=atoi(argv[4]), nprint=atoi(argv[8]), maxnn=atoi(argv[13]);
	bool kpts=atoi(argv[2]), pbc=atoi(argv[6]), v=atoi(argv[9]);
	int norbs=atoi(argv[12]), maxsteep=atoi(argv[16]);
	double rv=atof(argv[10]), rc=atof(argv[11]);
	double tol=atof(argv[15]), m=atof(argv[14]), dt=atof(argv[5]), T=atof(argv[7]), nu=atof(argv[17]), h=atof(argv[18]);
	std::vector<double> TBparam(6);
	TBparam[0]=atof(argv[19]);		// E_s
	TBparam[1]=atof(argv[20]);		// E_p
	TBparam[2]=atof(argv[21]);		// V_ss_sigma
	TBparam[3]=atof(argv[22]);		//	V_sp_sigma
	TBparam[4]=atof(argv[23]);		// V_pp_sigma
	TBparam[5]=atof(argv[24]);		// V_pp_pi
	std::vector<double> lats(3), posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz, &lats, pbc);
	int n=posx.size();
	// Determine maximum number of nearest neighbours
	if(n<maxnn){maxnn=n;}
	// Matrix neighbour list
	Eigen::MatrixXi inear(n,maxnn);
	std::vector<int> nnear(n);
	// Calculate distances
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	// Calculate distances
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n);
	Eigen::MatrixXd eigvects(norbs*n,norbs*n);
	// Geometry optimisation routine
	GeomOpt(norbs,rc,rv,m,dt,nmd,&posx,&posy,&posz,&refposx,&refposy,&refposz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr,&lats,pbc,T,nu,h,v,nprint,&TBparam,tol,maxsteep);

	return 0;
}

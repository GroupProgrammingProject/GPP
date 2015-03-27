#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "include/readinxyz.h"
#include "include/vectorfunctions.h"
#include "include/hamiltonian.h"
#include "include/functions.h"
#include "include/kpointsfunctions.h"
#include "include/geometryinfo.h"
#include "include/MolDyn.h"
#include "include/phonons.h"

int main(int argc, char* argv[]){
  	// Read in parameter values passed by the run script with MODE="energy"
	// xyzfilepath=argv[1], kptsfilepath=argv[3]
  //RUNCOMMAND="$0MAIN $1XYZ_FILE_PATH $2KPTS $3KPTS_FILE_PATH $4KSYMM $5NUM_STEPS $6DT $7PBC $8T $9FRAME_RATE $10VERBOSE $11RV $12RC $13NUM_ORBS $14MAX_NEIGHBOURS $15MASS $16TOL $17MAX_STEEP $18THERM_RATE $19H $20P1 $21P2 $22P3 $23P4 $24P5 $25P6"
	// From now on is the good version (input)
   int nmd=atoi(argv[5]), nprint=atoi(argv[9]), maxnn=atoi(argv[14]);
   bool kpts=atoi(argv[2]), pbc=atoi(argv[7]), v=atoi(argv[10]), ksymm=atoi(argv[4]);
	int norbs=atoi(argv[13]), maxsteep=atoi(argv[17]);
	double rv=atof(argv[11]), rc=atof(argv[12]);
	double tol=atof(argv[16]), m=atof(argv[15]), dt=atof(argv[6]), T=atof(argv[8]), nu=atof(argv[18]), h=atof(argv[19]);
	std::vector<double> TBparam(6);
	TBparam[0]=atof(argv[20]);		// E_s
	TBparam[1]=atof(argv[21]);		// E_p
	TBparam[2]=atof(argv[22]);		// V_ss_sigma
	TBparam[3]=atof(argv[23]);		//	V_sp_sigma
	TBparam[4]=atof(argv[24]);		// V_pp_sigma
	TBparam[5]=atof(argv[25]);		// V_pp_pi
	std::vector<double> lats(3), posx, posy, posz, vxin, vyin, vzin; 
	std::vector<bool> velspec;
	ReadInXYZ (argv[1], &posx, &posy, &posz, &vxin, &vyin, &vzin, &lats, pbc, &velspec);
	int i, n=posx.size();
	std::vector<double> fx(n), fy(n), fz(n), fmag(n);
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
//	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n);
	std::vector<double> refposx(n), refposy(n), refposz(n);
/*	velocity(m,&vx,&vy,&vz,T);
	//If input velocities are specified, initialise them for given atoms
	for(int i=0; i<n; i++){
		if(velspec.at(i)==1){
			vx.at(i)=vxin.at(i);
			vy.at(i)=vxin.at(i);
			vz.at(i)=vxin.at(i);
		}
	}
	//	std::vector<std::vector<double> > 
*/
	std::vector<std::pair<std::vector<double>,double> > kpoints; 
	Eigen::MatrixXd eigvects(norbs*n,norbs*n);                //real matrix
	if (kpts == 1) {
	  readinkpoints(argv[3],&kpoints,ksymm);
	}

/*
//Steepest descent relaxation routine
	double fmax,ebs,gam=1;
	int count=0,nmax=1e6; 
	std::cout << "Relaxing input structure..." << std::endl;
	do{
		forces(n,norbs,rc,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fx,&fy,&fz,&TBparam);
		for(i=0; i<n; i++){
			posx.at(i)=posx.at(i)+gam*h*fx.at(i);
			posy.at(i)=posy.at(i)+gam*h*fy.at(i);
			posz.at(i)=posz.at(i)+gam*h*fz.at(i);
			fmag.at(i)=sqrt(fx.at(i)*fx.at(i)+fy.at(i)*fy.at(i)+fz.at(i)*fz.at(i));
		}
		fmax=*std::max_element(fmag.begin(),fmag.end()); //find largest elemen of fmag
		GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
		NearestNeighbours(&inear,&nnear,&modr,rv);
		ebs=Hamiltonian(n,norbs,&TBparam,&modr,&rx,&ry,&rz,&eigvects,v);
		count=count+1;
	}while(fmax>1e-8 && count<nmax); //continue until desired accuracy reached, or we've reached nmax steps

	FILE *file_forces=fopen("relaxed_forces.xyz","w");
	for(i=0; i<n; i++){
		fprintf(file_forces,"%.10f %.10f %.10f\n", fx.at(i), fy.at(i), fz.at(i));
	}
	fclose(file_forces);

	FILE *file_rel=fopen("relaxed.xyz","w");
	fprintf(file_rel,"%d\nC%d\t%.3f \t %.3f \t %.3f \n",n,n,lats[0],lats[1],lats[2]);
	for(i=0; i<n; i++){
		fprintf(file_rel,"6 %.10f %.10f %.10f\n", posx.at(i), posy.at(i), posz.at(i));
	}
	fclose(file_rel);
*/


// Geometry optimisation routine
	GeomOpt(norbs,rc,rv,m,dt,nmd,&posx,&posy,&posz,&refposx,&refposy,&refposz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr,&lats,pbc,T,nu,h,v,nprint,&TBparam,tol,maxsteep,kpts,&kpoints);

	//Calculate normal modes
	normalmodes(n,norbs,rc,m,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&TBparam);

	return 0;
}

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "include/readinxyz.h"
#include "include/readinxyzv.h"
#include "include/vectorfunctions.h"
#include "include/hamiltonian.h"
#include "include/functions.h"
#include "include/geometryinfo.h"
#include "include/MolDyn.h"
#include <ctime>

int main(int argc, char* argv[]){
	// Read in the parameter values to run MD simulation
	// The order in which they are passed is important, determined in run script
	//RUNCOMMAND="$1MAIN $2XYZ_FILE_PATH $3NUM_STEPS $4DT $5PBC $6ENSEMBLE $7THERM_RATE $8T $9FRAME_RATE $10VERBOSE $11RV $12RC $13NUM_ORBS $14MAX_NEIGHBOURS $15MASS $16NU $17P1 $18P2 $19P3 $20P4 $21P5 $22P6"

	int nmd=atoi(argv[2]), nprint=atoi(argv[8]), norbs=atoi(argv[12]), maxnn=atoi(argv[13]);
	double dt=atof(argv[3]), T=atof(argv[7]), rv=atof(argv[10]),rc=atof(argv[11]), nu=atof(argv[6]),m=atof(argv[14]);
	bool v=atoi(argv[9]), pbc=atoi(argv[4]), ander=atoi(argv[5]);
	std::vector<double> TBparam(6); //vector to hold TB parameters (Xu et al)
	TBparam[0]=atof(argv[16]);		// E_s
	TBparam[1]=atof(argv[17]);		// E_p
	TBparam[2]=atof(argv[18]);		// V_ss_sigma
	TBparam[3]=atof(argv[19]);		//	V_sp_sigma
	TBparam[4]=atof(argv[20]);		// V_pp_sigma
	TBparam[5]=atof(argv[21]);		// V_pp_pi
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool renn;
	int i,j;
	std::vector<double> lats(3);
	// Read in types, 
	std::vector<double> posx, posy, posz, vxin, vyin, vzin;
	std::vector<bool> velspec;
//	ReadInXYZ (argv[1],&posx, &posy, &posz, &lats, pbc);
	ReadInXYZV (argv[1],&posx, &posy, &posz, &vxin, &vyin, &vzin, &lats, pbc, &velspec);
	// Number of atoms, number of orbitals, and number of MD steps
	int n=posx.size();
	// Velocities, reference postions, and vector neighbour list
	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n), fx(n), fy(n), fz(n);
	std::vector<int> nnear(n);
	// Determine maximum number of nearest neighbours
	if(n<maxnn){maxnn=n;}
	// Matrix neighbour list
	Eigen::MatrixXi inear(n,maxnn);
	// Calculate distances
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	// Timestep, initial temperature, atomic mass, cut off and Verlet radii
	double Tf,tmd,kb=1./11603;
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create empty arrays needed for MD

	//Matrix for eigenvectors
	Eigen::MatrixXd eigvects(norbs*n,norbs*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Calculation of nearest neighbours:
	NearestNeighbours(&inear,&nnear,&modr,rv);
	// Calculation of initial velocities:
	velocity(m,&vx,&vy,&vz,T);
	//If input velocities are given, initialise them for relevant atoms
	for(i=0; i<vx.size(); i++){
		if(velspec.at(i)==1){
			vx.at(i)=vxin.at(i);
			vy.at(i)=vyin.at(i);
			vz.at(i)=vzin.at(i);
		}
	}
	// Initialisation of reference positions:
	for(i=0;i<n;i++){
	  refposx.at(i)=posx.at(i);
	  refposy.at(i)=posy.at(i);
	  refposz.at(i)=posz.at(i);
	}
	FILE *file=fopen("movie.txt","w");
	fprintf(file,"%d\nC%d\n",n,n);
	for(i=0;i<n;i++){
	  fprintf(file,"6  %f %f %f \n",posx.at(i),posy.at(i),posz.at(i));
	}
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,norbs,&TBparam,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);	
	ekin=3*(n-1)*kb*T/2;
	etot=ebs+erep+ekin;
	FILE *en=fopen("energy.txt","w");
	FILE *f=fopen("forces.txt","w");
	fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",0.0,T,ekin,ebs,erep,etot);
	double xi1=0,xi2=0,vxi1=0,vxi2=0,q1=1,q2=1;
	// MD cycle
	for(i=0; i<nmd; i++){
	  Tf=verlet(norbs,rc,rv,m,dt,&posx,&posy,&posz,&refposx,&refposy,&refposz,&vx,&vy,&vz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr,ebs,&lats,pbc,T,nu,ander,&TBparam);
	  ekin=3*(n-1)*kb*Tf/2;
	  for(int k=0; k<n; k++){
			fprintf(f,"%f\t%f\t%f\t\n",fx.at(k),fy.at(k),fz.at(k));
	  }
	  if(i%nprint==0){
	    //H_MD and eigvects have now also been populated
	    erep=Erep(&modr);
	    etot=ebs+erep+ekin;
	    tmd=i*dt;
	    fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",tmd,Tf,ekin,ebs,erep,etot);
	    fprintf(file,"%d\nC%d molecule\n",n,n);
	    for(j=0;j<n;j++){
	      fprintf(file,"6  %f %f %f\n",posx.at(j),posy.at(j),posz.at(j));
	    }
	  } 
	  if(i*10%nmd==0){
	    std::cout << i*10/nmd << "0 % completed" << std::endl;
	  }
	}
	//Output final positions and velocities to a .xyz file for future simulations
	FILE *rel=fopen("final.xyz","w");
	fprintf(rel,"%d\nC%d molecule\n",n,n);
	for(j=0;j<n;j++){
		fprintf(rel,"6  %f %f %f %f %f %f\n",posx.at(j),posy.at(j),posz.at(j),vx.at(j),vy.at(j),vz.at(j));
	}
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,norbs,&TBparam,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
	etot=ebs+erep+ekin;

	std::cout << "Ebs = " << ebs << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	std::cout << "Etot = " << etot << std::endl;

return 0;
}

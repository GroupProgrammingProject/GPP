#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "include/readinxyz.h"
#include "include/vectorfunctions.h"
#include "include/hamiltonian.h"
#include "include/functions.h"
#include "include/geometryinfo.h"
#include "include/MolDyn.h"
#include <ctime>

int main(int argc, char* argv[]){
	// Read in the parameter values to run MD simulation
	// The order in which they are passed is important, determined in run script
	int nmd=atoi(argv[2]), nprint=atoi(argv[7]), norbs=atoi(argv[11]), maxnn=atoi(argv[12]);
	double dt=atof(argv[3]), T=atof(argv[6]), rv=atof(argv[9]),rc=atof(argv[10]);
	bool v=atoi(argv[8]), pbc=atoi(argv[4]), ander=atoi(argv[5]);
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool renn;
	int i,j;
	std::vector<double> lats(3);
	// Read in types, 
	std::vector<double> posx, posy, posz;
	std::vector<bool> velspec;
	ReadInXYZ (argv[1],&posx, &posy, &posz, &lats, pbc);
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
	double Tf,m=12*1.0365e2,tmd,kb=1./11603;
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create empty arrays needed for MD
//	std::cout << norbs << "  n=" << n << std::endl;

	std::vector<double> TBparam(6);
	// In final version, TB params will be read in from input file
	TBparam[0]=-2.99;		// E_s
	TBparam[1]=3.71;		// E_p
	TBparam[2]=-5;			// V_ss_sigma
	TBparam[3]=4.7;		//	V_sp_sigma
	TBparam[4]=5.5;		// V_pp_sigma
	TBparam[5]=-1.55;		// V_pp_pi

	//Matrix for eigenvectors
	Eigen::MatrixXd eigvects(norbs*n,norbs*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Calculation of nearest neighbours:
	NearestNeighbours(&inear,&nnear,&modr,rv);
	// Calculation of initial velocities:
	velocity(m,&vx,&vy,&vz,T);
	// Initialisation of reference positions:
	for(i=0;i<n;i++){
	  refposx.at(i)=posx.at(i);
	  refposy.at(i)=posy.at(i);
	  refposz.at(i)=posz.at(i);
	}
	FILE *file=fopen("movie.txt","w");
	fprintf(file,"%d\nC3 molecule\n",n);
	for(i=0;i<n;i++){
	  fprintf(file,"6  %f %f %f\n",posx.at(i),posy.at(i),posz.at(i));
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
	  forces(n,norbs,rc,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fx,&fy,&fz,&TBparam);
	  Tf=verlet(norbs,rc,rv,m,dt,&posx,&posy,&posz,&refposx,&refposy,&refposz,&vx,&vy,&vz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr,ebs,&lats,pbc,&TBparam);
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
	    fprintf(file,"%d\nC3 molecule\n",n);
	    for(j=0;j<n;j++){
	      fprintf(file,"6  %f %f %f\n",posx.at(j),posy.at(j),posz.at(j));
	    }
	  } 
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

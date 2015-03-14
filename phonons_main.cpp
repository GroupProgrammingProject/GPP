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
#include "include/phonons.h"

int main(int argc, char* argv[]){
	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0, renn;
	int i,j;
	std::vector<double> lats(3), posx, posy, posz;
	bool pbc = 1;
	ReadInXYZ (argv[1],&posx, &posy, &posz, &lats, pbc);
//	scramble(&posx,&posy,&posz);
	// Number of atoms, number of orbitals, max no. of neares neighbours and number of MD steps
	int n=posx.size(),norbs=4,nmd=1,nprint=1;
	// Determine maximum number of nearest neighbours
	int maxnn=100;
	if(n<maxnn){
	  maxnn=n;
	}
	// Velocities, reference postions, and vector neighbour list
	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n),fx(n),fy(n),fz(n),fmag(n);
	// Matrix neighbour list
	std::vector<int> nnear(n);
	Eigen::MatrixXi inear(n,maxnn);
	// Timestep, initial temperature, atomic mass, cut off and Verlet radii
	double dt=1,T=500,Tf,m=12*1.0365e2,rc=2.6,rv=3,tmd,kb=1./11603,fmax,h=0.001,gam=1;
	int nmax=10000,count=0;
	// Calculate distances
	Eigen::MatrixXd modr(n,n), rx(n,n), ry(n,n), rz(n,n);
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	//TB parameters
	std::vector<double> TBparam(6);
	// In final version, TB params will be read in from input file
	TBparam[0]=-2.99;		// E_s
	TBparam[1]=3.71;		// E_p
	TBparam[2]=-5;			// V_ss_sigma
	TBparam[3]=4.7;		//	V_sp_sigma
	TBparam[4]=5.5;		// V_pp_sigma
	TBparam[5]=-1.55;		// V_pp_pi

// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Calculation of nearest neighbours:
	NearestNeighbours(&inear,&nnear,&modr,rv);
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,norbs,&TBparam,&modr,&rx,&ry,&rz,&eigvects,v);

	//Steepest descent relaxation routine
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

	FILE *file_rel=fopen("relaxed.xyz","w");
	fprintf(file_rel,"%d\nC%d\t%.3f \t %.3f \t %.3f \n",n,n,lats[0],lats[1],lats[2]);
	for(i=0; i<n; i++){
		fprintf(file_rel,"6 %.10f %.10f %.10f\n", posx.at(i), posy.at(i), posz.at(i));
	}
	fclose(file_rel);

	//Calculate normal modes
	std::cout << "Calculating normal modes..." << std::endl;
	std::vector<double> eigfreq(3*n);
	normalmodes(n,norbs,rc,m,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fx,&fy,&fz,&eigfreq,&TBparam);
	std::sort (eigfreq.begin(), eigfreq.end());
	std::reverse (eigfreq.begin(), eigfreq.end());
	FILE *file_freq=fopen("vibrational_frequencies.dat","w");
//std::cout << "Real eigenvalues expressed as wavectors in cm-1" << std::endl;
//for(i=0;i<3*n;i++){	std::cout << eigfreq[i]/n << std::endl; }
	for(i=0;i<3*n;i++){
		fprintf(file_freq,"%f\n", eigfreq[i]/n);
	}
	fclose(file_freq);

return 0;
}

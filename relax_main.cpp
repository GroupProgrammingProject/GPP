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

	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0, renn;
	int i,j;
	std::vector<double> lats(3);
	// Read in types, 
	std::vector<double> posx, posy, posz;
	bool pbc = 1;
	ReadInXYZ (argv[1],&posx, &posy, &posz, &lats, pbc);
//	scramble(&posx,&posy,&posz);
	// Number of atoms, number of orbitals, and number of MD steps
	int n=posx.size(),norbs=4,nmd=1,nprint=1;
	// Velocities, reference postions, and vector neighbour list
	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n),fx(n),fy(n),fz(n),fmag(n);
	std::vector<int> nnear(n);
	// Determine maximum number of nearest neighbours
	int maxnn=100;
	if(n<maxnn){
	  maxnn=n;
	}
	// Matrix neighbour list
	Eigen::MatrixXi inear(n,maxnn);
	// Calculate distances
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	// Timestep, initial temperature, atomic mass, cut off and Verlet radii
	double dt=1,T=500,Tf,m=12*1.0365e2,rc=2.6,rv=3,tmd,kb=1./11603,fmax,h=0.001,gam=1;
	int nmax=10000,count=0;
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Calculation of nearest neighbours:
	NearestNeighbours(&inear,&nnear,&modr,rv);
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	FILE *file=fopen("movie.txt","w");
	FILE *file2=fopen("forces.txt","w");
	fprintf(file,"%d\nC12 \t %.3f \t %.3f \t %.3f \n",n,lats[0],lats[1],lats[2]);
	for(i=0; i<n; i++){
		fprintf(file,"6 %.10f %.10f %.10f\n", posx.at(i), posy.at(i), posz.at(i));
	}
	//steepest descent relaxation routine
	do{
		if(count*100%nmax==0){std::cout << count*100/nmax << "% complete" << std::endl;}
		forces(n,norbs,rc,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fx,&fy,&fz);
		fprintf(file,"%d\nC12 \t %.3f \t %.3f \t %.3f \n",n,lats[0],lats[1],lats[2]);
		for(i=0; i<n; i++){
			posx.at(i)=posx.at(i)+gam*h*fx.at(i);
			posy.at(i)=posy.at(i)+gam*h*fy.at(i);
			posz.at(i)=posz.at(i)+gam*h*fz.at(i);
			fmag.at(i)=sqrt(fx.at(i)*fx.at(i)+fy.at(i)*fy.at(i)+fz.at(i)*fz.at(i));
			fprintf(file,"6 %.10f %.10f %.10f\n", posx.at(i), posy.at(i), posz.at(i));
			fprintf(file2,"6 %i %.10f %.10f %.10f %.10f \n", count, fx.at(i), fy.at(i), fz.at(i), fmag.at(i));
		}
		fmax=*std::max_element(fmag.begin(),fmag.end()); //find largest elemen of fmag
		GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
		NearestNeighbours(&inear,&nnear,&modr,rv);
		ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
		count=count+1;
	}while(fmax>1e-8 && count<nmax); //continue until desired accuracy reached, or we've reached nmax steps
	FILE *file_rel=fopen("relax.txt","w");
	fprintf(file_rel,"%d\nC12 \t %.3f \t %.3f \t %.3f \n",n,lats[0],lats[1],lats[2]);
	for(i=0; i<n; i++){
		fprintf(file_rel,"6 %.10f %.10f %.10f\n", posx.at(i), posy.at(i), posz.at(i));
	}

return 0;
}

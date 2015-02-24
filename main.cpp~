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

int main(int argc, char* argv[]){
	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0, renn;
	int i,j;
	// Read in types, 
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1],&posx, &posy, &posz);
	// Number of atoms, number of orbitals, and number of MD steps
	int n=posx.size(),norbs=4,nmd=1000,nprint=10;
	// Velocities, reference postions, and vector neighbour list
	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n);
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
	GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz);
	// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Timestep, initial temperature, atomic mass, cut off and Verlet radii
	double dt=1,T=1000,Tf,m=12*1.0365e2,rc=2.6,rv=3,tmd,kb=1./11603;
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
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);	
	ekin=3*n*kb*T/2;
	etot=ebs+erep+ekin;

	FILE *en=fopen("energy.txt","w");
	fprintf(en,"%f\t%f\n",0.0,etot);

	
	// MD cycle
	for(i=0;i<nmd;i++){
	  Tf=verlet(norbs,rc,rv,m,dt,&posx,&posy,&posz,&refposx,&refposy,&refposz,&vx,&vy,&vz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr);
	  ekin=3*n*kb*Tf/2;
	  
	  if(i%nprint==0){
	    // Starting TB module: calculating energies
	    ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	    //H_MD and eigvects have now also been populated
	    erep=Erep(&modr);
	    etot=ebs+erep+ekin;
	    tmd=i*dt;
	    fprintf(en,"%f\t%f\n",tmd,etot);
	    fprintf(file,"%d\nC3 molecule\n",n);
	    for(j=0;j<n;j++){
	      fprintf(file,"6  %f %f %f\n",posx.at(j),posy.at(j),posz.at(j));
	    }
	  } 
	}

	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
	etot=ebs+erep;

	std::cout << "Ebs = " << ebs << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	std::cout << "Etot = " << etot << std::endl;

return 0;
}

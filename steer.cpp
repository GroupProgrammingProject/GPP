// File that assembles all the different TB and MD modules
#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "readinxyz.h"
#include "vectorfunctions.h"
#include "hamiltonian.h"
#include "functions.h"
#include "MolDyn.h"

int main(int argc, char* argv[]){

	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	// Read in types, 
	std::vector<int> type;
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &type,&posx, &posy, &posz);
	// Number of atoms
	int n=posx.size();
	// Create empty arrays needed for MD
	std::vector<double> eigvects(16*n*n);
	// Energies from TB model
	double ebs,erep,etot;

	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&posx,&posy,&posz,&eigvects);
	//H_MD and eigvects have now also been populated
	erep=Erep(&posx,&posy,&posz);
// Determining erep works, however, we need to check if we're passing pointers or arrays to Erep()
	etot=ebs+erep;

	// Start MD routine: steepest descent
	// First call near_neigh()
	int i,j,nnmax=10;	
	std::vector<double> mass(n);
	std::vector<int> nnear(n),inear(n*nnmax);		// Mass in kg?
	double m=12;
	for(i=0;i<n;i++){mass.at(i)=m;}
	double rc=2.6;
	double sx=10,sy=10,sz=10;
	near_neigh(n,&posx,&posy,&posz,rc,&nnear,&inear,sx,sy,sz);
//test near_neigh
//	for(i=0;i<n;i++){std::cout << "inear of atom " << i << " is " << inear.at(i) << std::endl;}

	// Start forces convergence loop
	int norbs=4;
	std::vector<double> fx(n),fy(n),fz(n);
	double h=0.01;
	
	//	FILE *start = fopen("start.txt", "w");
	//for(int i=0;i<n;i++){
	//fprintf(start,"%f\t%f\t%f\n",posx.at(i),posy.at(i),posz.at(i));
	//}

	//for(int i=0;i<10000;i++){
	  forces(n,norbs,&posx,&posy,&posz,&eigvects,rc,&nnear,&inear,&fx,&fy,&fz);
	  /*	  for(int j=0;j<n;j++){
	    posx.at(j)=posx.at(j)+h*fx.at(j);
	    posy.at(j)=posy.at(j)+h*fy.at(j);
	    posz.at(j)=posz.at(j)+h*fz.at(j);
	  }
	}

	/*	FILE *end = fopen("end.txt", "w");
	for(int i=0;i<n;i++){
	  fprintf(end,"%f\t%f\t%f\n",posx.at(i),posy.at(i),posz.at(i));
	}	
       
	/*     	for(int i=0;i<n;i++){
	  std::cout << "fx(" << i << ")=" << fx.at(i) << std::endl;
	  std::cout << "fy(" << i << ")=" << fy.at(i) << std::endl;
	  std::cout << "fz(" << i << ")=" << fz.at(i) << std::endl;
	  std::cout << std::endl;
	}
	*/
//	change positions of atoms

	//Plot positions!
	

return 0;
}

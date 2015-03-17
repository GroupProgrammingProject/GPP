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
	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0, renn;
	int i,j;
	std::vector<double> lats(3);
	// Read in types, 
	std::vector<double> posx, posy, posz, vxin, vyin, vzin;
	bool pbc = 1, ander=0;
	std::vector<bool> velspec;
	ReadInXYZV (argv[1],&posx, &posy, &posz, &vxin, &vyin, &vzin, &lats, pbc,&velspec);
	// Number of atoms, number of orbitals, and number of MD steps
	int n=posx.size(),norbs=4,nmd=1000,nprint=1;
	// Velocities, reference postions, and vector neighbour list
	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n), fx(n), fy(n), fz(n);
	std::vector<int> nnear(n);
	// Determine maximum number of nearest neighbours
	int maxnn=100;
	if(n<maxnn){maxnn=n;}
	// Matrix neighbour list
	Eigen::MatrixXi inear(n,maxnn);
	// Calculate distances
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	// Timestep, initial temperature, atomic mass, cut off and Verlet radii
	double dt=1,T=100,Tf,m=12*1.0365e2,rc=2.6,rv=3,tmd,kb=1./11603,nu=0.5;
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create empty arrays needed for MD
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
	fprintf(file,"%d\nC3 molecule\n",n);
	for(i=0;i<n;i++){
	  fprintf(file,"6  %f %f %f %f %f %f\n",posx.at(i),posy.at(i),posz.at(i),vx.at(i),vy.at(i),vz.at(i));
	}
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);	
	ekin=3*(n-1)*kb*T/2;
	etot=ebs+erep+ekin;
	FILE *en=fopen("energy.txt","w");
	FILE *f=fopen("forces.txt","w");
	fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",0.0,T,ekin,ebs,erep,etot);
	double xi1=0,xi2=0,vxi1=0,vxi2=0,q1=1,q2=1;
	// MD cycle
	for(i=1;i<nmd+1;i++){
		if(i*10%nmd==0){
			std::cout << i*10/nmd << "0% completed" << std::endl;}
	  forces(n,norbs,rc,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fx,&fy,&fz);
	  Tf=verlet(norbs,rc,rv,m,dt,&posx,&posy,&posz,&refposx,&refposy,&refposz,&vx,&vy,&vz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr,ebs,&lats,pbc,T,nu,ander);
	  //canonical ensemble function
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
			if(pbc==1){fprintf(file,"%d\nC%d.xyz %f %f %f\n",n,n,lats.at(0),lats.at(1),lats.at(2));}
			else{fprintf(file,"%d\nC%d.xyz\n",n,n);}
	    	for(j=0;j<n;j++){
	      	fprintf(file,"6  %f %f %f %f %f %f\n",posx.at(j),posy.at(j),posz.at(j),vx.at(j), vy.at(j), vz.at(j));
			}
	  } 
	}
	FILE *fin=fopen("final.xyz","w");
	if(pbc==1){fprintf(fin,"%d\nC%d.xyz %f %f %f\n",n,n,lats.at(0),lats.at(1),lats.at(2));}
	else{fprintf(fin,"%d\nC%d.xyz\n",n,n);}
	//Output final positions and velocities to a .xyz file for future simulations
	for(j=0;j<n;j++){
		fprintf(fin,"6  %f %f %f %f %f %f\n",posx.at(j),posy.at(j),posz.at(j),vx.at(j), vy.at(j), vz.at(j));
	}
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
	etot=ebs+erep+ekin;

	std::cout << "Ebs = " << ebs << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	std::cout << "Etot = " << etot << std::endl;
	std::cout << lats.at(0) << " " << lats.at(1) << " " << lats.at(2) << std::endl;

return 0;
}

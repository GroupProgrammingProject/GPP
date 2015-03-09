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

void normalmodes(int n,int norbs,double rc,double m,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* eigvects, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz, std::vector<double>* eigfreq);
/*Inputs, in order: #atoms; #orbitals; cut-off radius; mass of atoms;atom vector distances; modulus of vector distances; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; forces vectors,eigenmodes*/

void normalmodes(int n,int norbs,double rc,double m,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* eigvects, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz, std::vector<double>* eigfreq){

	bool v=0;
	int i,j,l,J;						//J is index of atom to be displaced
	int beta;							//direction of displacement (0=x, 1=y, 2=z)
	double diff, abs_diff=1e-5;	//distance of displacement
	double ebs;							//need for calling Hamiltonian
	double h_bar=6.582119e-16;		// eV*s
	double c=29979245800;			// cm/s
	Eigen::MatrixXd modr_0(n,n);
	Eigen::MatrixXd rx_0(n,n);
	Eigen::MatrixXd ry_0(n,n);
	Eigen::MatrixXd rz_0(n,n);
	Eigen::MatrixXd fr(3*n,3*n), fl(3*n,3*n);
	Eigen::MatrixXd dynamicmat(3*n,3*n);

 	ebs=Hamiltonian(n,modr,rx,ry,rz,eigvects,v);
	forces(n,norbs,rc,rx,ry,rz,modr,eigvects,nnear,inear,fx,fy,fz);//recalculate forces
	
std::cout << "forces of (hopefully) equilibrium structure:" << std::endl;
for(i=0;i<n;i++){
	std::cout << "fx = " << fx->at(i) << std::endl;
	std::cout << "fy = " << fy->at(i) << std::endl;
	std::cout << "fz = " << fz->at(i) << std::endl << std::endl;
}



	//store original distances
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			rx_0(i,j)=(*rx)(i,j);
			ry_0(i,j)=(*ry)(i,j);
			rz_0(i,j)=(*rz)(i,j);
			modr_0(i,j)=(*modr)(i,j);
		}
	}
	//set dynamicmat, fr and fl to zero
	for(i=0;i<3*n;i++){
		for(j=0;j<3*n;j++){
			dynamicmat(i,j)=0;
			fr(i,j)=0;
			fl(i,j)=0;
		}
	}	

	// Start loop for displacing atoms in positive (l=0) and negative (l=1) directions
	for(l=0;l<2;l++){
		if(l==0) { diff=-abs_diff; }			//fill up fr (moving atom J in beta direction by diff)
		else if(l==1) { diff=abs_diff; } 		//fill up fl (move atom by -diff)
		//Start loop over all n atoms
		for(J=0;J<n;J++){
			//Start loop over all three directions
			for(beta=0;beta<3;beta++){
				//revert all distances to original	
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						(*rx)(i,j)=rx_0(i,j);
						(*ry)(i,j)=ry_0(i,j);
						(*rz)(i,j)=rz_0(i,j);
						(*modr)(i,j)=modr_0(i,j);
					}
				}

/*		std::cout << "original geometry" << std::endl;
		std::cout << "rx" << std:: endl << rx << std::endl;
		std::cout << "ry" << std:: endl << ry << std::endl;
		std::cout << "rz" << std:: endl << rz << std::endl;	
		std::cout << "modr" << std:: endl << modr << std::endl;
*/
				//Update all distances (keep its own to 0)
				for(i=0;i<n;i++){
					if(i!=J){
						if(beta==0){ 
							(*rx)(i,J)=(*rx)(i,J)+diff;
							(*rx)(J,i)=-(*rx)(i,J);
						}
						else if (beta==1){ 
							(*ry)(i,J)=(*ry)(i,J)+diff;
							(*ry)(J,i)=-(*ry)(i,J);
						}
						else if (beta==2){ 
							(*rz)(i,J)=(*rz)(i,J)+diff; 
							(*rz)(J,i)=-(*rz)(i,J);
						}
					}
				}	
				for(i=0;i<n;i++){
					if(i!=J){
						(*modr)(i,J)=sqrt( (*rx)(i,J) * (*rx)(i,J) + (*ry)(i,J) * (*ry)(i,J) + (*rz)(i,J) * (*rz)(i,J) );
					}
				}
/*
std::cout << "Moved atom " << J << " in " << beta << " direction" << std::endl;
std::cout << "rx" << std:: endl << rx << std::endl;
std::cout << "ry" << std:: endl << ry << std::endl;
std::cout << "rz" << std:: endl << rz << std::endl;	
std::cout << "modr" << std:: endl << modr << std::endl;
*/
				// Update eigenvectors and forces for new distances
		 		ebs=Hamiltonian(n,modr,rx,ry,rz,eigvects,v);
				forces(n,norbs,rc,rx,ry,rz,modr,eigvects,nnear,inear,fx,fy,fz);//recalculate forces
/*	
std::cout << "new forces:" << std::endl;
for(i=0;i<n;i++){
	std::cout << "fx = " << fx->at(i) << std::endl;
	std::cout << "fy = " << fy->at(i) << std::endl;
	std::cout << "fz = " << fz->at(i) << std::endl << std::endl;
}
*/
				for(i=0;i<n;i++){
					if(l==0){
						fr(3*i,3*J+beta)=fx->at(i);
						fr(3*i+1,3*J+beta)=fy->at(i);
						fr(3*i+2,3*J+beta)=fz->at(i);
					}
					else if (l==1){
						fl(3*i,3*J+beta)=fx->at(i);
						fl(3*i+1,3*J+beta)=fy->at(i);
						fl(3*i+2,3*J+beta)=fz->at(i);
					}
				}
			}	//end beta loop
		}	//end J loop
	}	//end l loop

/*
std::cout << "Final fr: " << std::endl;
std::cout << fr << std::endl << std::endl;
std::cout << "Final fl " << std::endl;
std::cout << fl << std::endl;
*/

	//Create dynamical matrix
	for(i=0;i<3*n;i++){
		for(j=0;j<3*n;j++){
			dynamicmat(i,j)=(fl(i,j)-fr(i,j))/(2*abs_diff*m);
		}
	}

//std::cout << "Dynamical matrix:" << std::endl;
//std::cout << dynamicmat << std::endl << std::endl;

	//Solve for eigenmodes and eigenfrequencies
	Eigen::EigenSolver<Eigen::MatrixXd> es(dynamicmat);
	//returns eigenfrequencies in fs-1

std::cout << "eigenvalues:" << std::endl;
std::cout << es.eigenvalues() << std::endl;

//std::cout << std::endl << "eigenvectors:" << std::endl;
//std::cout << es.eigenvectors() << std::endl;

	for(i=0;i<3*n;i++){ (*eigfreq).at(i) = 1e15*sqrt(real(es.eigenvalues()[i]))/c; }
	//convert to wavevectors in cm-1

}	//end eigenmodes()

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
	// Number of atoms, number of orbitals, and number of MD steps
	int n=posx.size(),norbs=4,nmd=1,nprint=1;
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
	// Timestep, initial temperature, atomic mass, cut off and Verlet radii
	double dt=1,T=500,Tf,m=12*1.0365e2,rc=2.6,rv=3,tmd,kb=1./11603;
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Calculation of nearest neighbours:
	NearestNeighbours(&inear,&nnear,&modr,rv);

	erep=Erep(&modr);	
/*	Erep has to be called, otherwise errors appear: several functions from "functions.h" are undefined"	*/

	//Calculate eigenmodes and eigenfrequencies
	std::vector<double> fx(n),fy(n),fz(n),eigfreq(3*n);
	normalmodes(n,norbs,rc,m,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fx,&fy,&fz,&eigfreq);

std::cout << "Real eigenvalues expressed as wavectors in cm-1" << std::endl;
for(i=0;i<3*n;i++){	std::cout << eigfreq[i] << std::endl; }

return 0;
}

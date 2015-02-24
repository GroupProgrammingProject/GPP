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
	std::vector<double> vx(n), vy(n), vz(n), refposx(n), refposy(n), refposz(n), fx(n), fy(n), fz(n), dposx(n), dposy(n), dposz(n);
	std::vector<int> nnear(n);
	// Determine maximum number of nearest neighbours
	int maxnn=100;
	if(n<maxnn){
	  maxnn=n;
	}
	// Matrix neighbour list
	Eigen::MatrixXi inear(n,maxnn);
	// Calculate distances
	Eigen::MatrixXd modr(n,n),drmodr(n,n),dlmodr(n,n);
	Eigen::MatrixXd rx(n,n),drrx(n,n),dlrx(n,n);
	Eigen::MatrixXd ry(n,n),drry(n,n),dlry(n,n);
	Eigen::MatrixXd rz(n,n),drrz(n,n),dlrz(n,n);
	GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz);
	// Create empty arrays needed for MD
	Eigen::MatrixXd eigvects(4*n,4*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Timestep, initial temperature, atomic mass, cut off and Verlet radii
	double dt=1,T=1000,Tf,m=12*1.0365e2,rc=2.6,rv=3,tmd,kb=1./11603,fbsx,fbsy,fbsz,frepx,frepy,frepz,h=0.1;
	// Calculation of nearest neighbours:
	NearestNeighbours(&inear,&nnear,&modr,rv);
	Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
      
	forces(n,norbs,rc,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fx,&fy,&fz);

	for(i=0;i<n;i++){
	  std::cout << "forces routine; i=" << i << std::endl;
	  std::cout << fx.at(i) << " " << fy.at(i) << " " << fz.at(i) << std::endl;
	  
	}
	
	for(i=0;i<n;i++){
	  dposx.at(i)=posx.at(i);
	  dposy.at(i)=posy.at(i);
	  dposz.at(i)=posz.at(i);
	}

	for(i=0;i<n;i++){
	  std::cout << "i=" << i << std::endl;
	  
	  dposx.at(i)=posx.at(i)+h;
	  GetAllDistances(&drmodr,&drrx,&drry,&drrz,&dposx,&posy,&posz);
	  dposx.at(i)=posx.at(i)-h;
	  GetAllDistances(&dlmodr,&dlrx,&dlry,&dlrz,&dposx,&posy,&posz);
	  fbsx=-(Hamiltonian(n,&drmodr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&dlmodr,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
	  frepx=-(Erep(&drmodr)-Erep(&dlmodr))/(2*h);

	  dposy.at(i)=posy.at(i)+h;
	  GetAllDistances(&drmodr,&drrx,&drry,&drrz,&posx,&dposy,&posz);
	  dposy.at(i)=posy.at(i)-h;
	  GetAllDistances(&dlmodr,&dlrx,&dlry,&dlrz,&posx,&dposy,&posz);
	  fbsy=-(Hamiltonian(n,&drmodr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&dlmodr,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
	  frepy=-(Erep(&drmodr)-Erep(&dlmodr))/(2*h);
	  
	  dposz.at(i)=posz.at(i)+h;
	  GetAllDistances(&drmodr,&drrx,&drry,&drrz,&posx,&posy,&dposz);
	  dposz.at(i)=posz.at(i)-h;
	  GetAllDistances(&dlmodr,&dlrx,&dlry,&dlrz,&posx,&posy,&dposz);
	  fbsz=-(Hamiltonian(n,&drmodr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&dlmodr,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
	  frepz=-(Erep(&drmodr)-Erep(&dlmodr))/(2*h);
	  
	  std::cout << Erep(&drmodr) << " " << Erep(&dlmodr) << std::endl;

	  for(j=0;j<n;j++){
	    dposx.at(j)=posx.at(j);
	    dposy.at(j)=posy.at(j);
	    dposz.at(j)=posz.at(j);
	  }

	  
	  std::cout << fbsx << " " << fbsy << " " << fbsz << std::endl;
	  
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

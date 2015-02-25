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
	int n=posx.size(),norbs=4,nmd=1000,nprint=1;
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
	double dt=1,T=500,Tf,m=12*1.0365e2,rc=2.6,rv=3,tmd,kb=1./11603;
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

	Eigen::MatrixXd modrl(n,n);
	Eigen::MatrixXd dlrx(n,n);
	Eigen::MatrixXd dlry(n,n);
	Eigen::MatrixXd dlrz(n,n);
	Eigen::MatrixXd modrr(n,n);
	Eigen::MatrixXd drrx(n,n);
	Eigen::MatrixXd drry(n,n);
	Eigen::MatrixXd drrz(n,n);
	std::vector<double> drposx(n),drposy(n),drposz(n),f_fdx(n),f_fdy(n),f_fdz(n),fxa(n),fya(n),fza(n),df1x(n),df1y(n),df1z(n),df2x(n),df2y(n),df2z(n),f_repx(n),f_repy(n),f_repz(n);
	std::vector<double> dlposx(n),dlposy(n),dlposz(n);
	double h=0.0001;
	for(int k=0; k<n; k++)
	{
		drposx.at(k)=posx.at(k);
		drposy.at(k)=posy.at(k);
		drposz.at(k)=posz.at(k);
		dlposx.at(k)=posx.at(k);
		dlposy.at(k)=posy.at(k);
		dlposz.at(k)=posz.at(k);
	}

	FILE *en=fopen("energy.txt","w");
	fprintf(en,"%f\t%f\n",0.0,etot);
	// MD cycle
	for(i=0;i<nmd;i++){
//		forces(n,norbs,rc,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fxa,&fya,&fza);
		Tf=verlet(norbs,rc,rv,m,dt,&posx,&posy,&posz,&refposx,&refposy,&refposz,&vx,&vy,&vz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr);
		GetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz);
		forces(n,norbs,rc,&rx,&ry,&rz,&modr,&eigvects,&nnear,&inear,&fxa,&fya,&fza);
		ekin=3*n*kb*Tf/2;
		if(i%nprint==0){ 
			std::cout << "i=" << i << std::endl << std::endl;
			for(int k=0; k<n; k++){
				for(int j=0; j<n; j++)
				{
					drposx.at(j)=posx.at(j);
					drposy.at(j)=posy.at(j);
					drposz.at(j)=posz.at(j);
					dlposx.at(j)=posx.at(j);
					dlposy.at(j)=posy.at(j);
					dlposz.at(j)=posz.at(j);
				}
				drposx.at(k)=posx.at(k)+h;
				dlposx.at(k)=posx.at(k)-h;
				GetAllDistances(&modrr,&drrx,&drry,&drrz,&drposx,&posy,&posz);
				GetAllDistances(&modrl,&dlrx,&dlry,&dlrz,&dlposx,&posy,&posz);
				f_fdx.at(k)=-(Hamiltonian(n,&modrr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&modrl,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
				f_repx.at(k)=-(Erep(&modrr)-Erep(&modrl))/(2*h);
//				df1x.at(k)=fabs(fxa.at(k)-f_fdx.at(k));
				
				drposy.at(k)=posy.at(k)+h;
				dlposy.at(k)=posy.at(k)-h;
				GetAllDistances(&modrr,&drrx,&drry,&drrz,&posx,&drposy,&posz);
				GetAllDistances(&modrl,&dlrx,&dlry,&dlrz,&posx,&dlposy,&posz);
				f_fdy.at(k)=-(Hamiltonian(n,&modrr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&modrl,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
				f_repy.at(k)=-(Erep(&modrr)-Erep(&modrl))/(2*h);
//				df1y.at(k)=fabs(fya.at(k)-f_fdy.at(k));
				
				drposz.at(k)=posz.at(k)+h;
				dlposz.at(k)=posz.at(k)-h;
				GetAllDistances(&modrr,&drrx,&drry,&drrz,&posx,&posy,&drposz);
				GetAllDistances(&modrl,&dlrx,&dlry,&dlrz,&posx,&posy,&dlposz);
				f_fdz.at(k)=-(Hamiltonian(n,&modrr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&modrl,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
				f_repz.at(k)=-(Erep(&modrr)-Erep(&modrl))/(2*h);
//				df1z.at(k)=fabs(fza.at(k)-f_fdz.at(k));
				std::cout << "Analytic forces:  " << fxa.at(k) << "  " << fya.at(k) << "  " << fza.at(k) << std::endl;
				std::cout << "Finite difference (BS):  " << f_fdx.at(k) << "  " << f_fdy.at(k) << "  " << f_fdz.at(k) << std::endl;
//				std::cout << "Analytic forces (R)" << fxa.at(k) << "  " << fya.at(k) << "  " << fza.at(k) << std::endl;
//				std::cout << "Finite difference (R):  " << f_repx.at(k) << "  " << f_repy.at(k) << "  " << f_repz.at(k) << std::endl;
//				std::cout << "Total FD force:  " << f_repx.at(k)+f_fdx.at(k) << "  " << f_repy.at(k)+f_fdy.at(k) << "  " << f_repz.at(k)+f_fdz.at(k) << std::endl;
			}
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
	etot=ebs+erep+ekin;

	std::cout << "Ebs = " << ebs << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	std::cout << "Ekin = " << ekin << std::endl;
	std::cout << "Etot = " << etot << std::endl;

return 0;
}

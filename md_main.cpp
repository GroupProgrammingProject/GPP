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
	//RUNCOMMAND="$0MAIN $1XYZ_FILE_PATH $2NUM_STEPS $3DT $4PBC $5ENSEMBLE $6THERM_RATE $7T $8FRAME_RATE $9VERBOSE $10RV $11RC $12NUM_ORBS $13MAX_NEIGHBOURS $14MASS $15P1 $16P2 $17P3 $18P4 $19P5 $20P6"

	//nmd= no. of MD steps, nprint=how often steps are printed to file, norbs=no. of atomic orbitals, maxnn=max no. of nearest neighbours that any atom can have.
	int nmd=atoi(argv[2]), nprint=atoi(argv[8]), norbs=atoi(argv[12]), maxnn=atoi(argv[13]);
	//dt=time step, T=temperature, rv=Verlet radius, rc=cut-off radius (from Xu et al), nu=frequency of collisions in canonical ensemble, m =mass.
	double dt=atof(argv[3]), T=atof(argv[7]), rv=atof(argv[10]), rc=atof(argv[11]), nu=atof(argv[6]), m=atof(argv[14]);
	//v=verbose (turn Hamiltonian printing on/off), pbc=switch for PBCs, ander=switch for thermof=dynamic ensemble.
	bool v=atoi(argv[9]), pbc=atoi(argv[4]), ander=atoi(argv[5]);
	std::vector<double> TBparam(6); //vector to hold TB parameters
	TBparam[0]=atof(argv[15]);		// E_s
	TBparam[1]=atof(argv[16]);		// E_p
	TBparam[2]=atof(argv[17]);		// V_ss_sigma
	TBparam[3]=atof(argv[18]);		//	V_sp_sigma
	TBparam[4]=atof(argv[19]);		// V_pp_sigma
	TBparam[5]=atof(argv[20]);		// V_pp_pi
	
	bool renn; //tells us whether or not we need to recalculate nearest neighbours at each MD step
	int i,j; //dummy variables for looping
	std::vector<double> lats(3);
	// Read in types, 
	std::vector<double> posx, posy, posz, vxin, vyin, vzin;
	std::vector<bool> velspec;
	//Read in Cartesian positions of atoms + lattice parameters and velocities (if applicable)
	ReadInXYZ (argv[1],&posx, &posy, &posz, &vxin, &vyin, &vzin, &lats, pbc, &velspec);
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
	double Tf,tmd,kb=1./11603,msvx,msvy,msvz;
	GetDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,&lats,rv,pbc);
	//Matrix for eigenvectors
	Eigen::MatrixXd eigvects(norbs*n,norbs*n);
	// Energies from TB model
	double ebs,erep,etot,ekin;
	// Calculation of nearest neighbours:
	NearestNeighbours(&inear,&nnear,&modr,rv);
	// Calculation of initial velocities
	velocity(m,&vx,&vy,&vz,T);
	//If input velocities are given, initialise them for relevant atoms
	for(i=0; i<n; i++){
		if(velspec.at(i)==1){
			vx.at(i)=vxin.at(i);
			vy.at(i)=vyin.at(i);
			vz.at(i)=vzin.at(i);
		}
		msvx=msvx+vx.at(i)*vx.at(i);
		msvy=msvy+vy.at(i)*vy.at(i);
		msvz=msvz+vz.at(i)*vz.at(i);
	}
	// Initialisation of reference positions
	for(i=0;i<n;i++){
	  refposx.at(i)=posx.at(i);
	  refposy.at(i)=posy.at(i);
	  refposz.at(i)=posz.at(i);
	}
	//Open file to store atom positions in MD steps
	FILE *mov=fopen("movie.txt","w");
	fprintf(mov,"%d\nC%d\n",n,n);
	for(i=0;i<n;i++){
	  fprintf(mov,"6  %f %f %f \n",posx.at(i),posy.at(i),posz.at(i));
	}
	// Starting TB	module: calculating energies
	ebs=Hamiltonian(n,norbs,&TBparam,&modr,&rx,&ry,&rz,&eigvects,v);
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);	
	ekin=0.5*m*(msvx+msvy+msvz);
	etot=ebs+erep+ekin;
	//Open files to store energies at each step
	FILE *en=fopen("energy.txt","w");
	fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",0.0,T,ekin,ebs,erep,etot);
	double xi1=0,xi2=0,vxi1=0,vxi2=0,q1=1,q2=1;
	// MD cycle
	for(i=1; i<nmd+1; i++){
	  Tf=verlet(norbs,rc,rv,m,dt,&posx,&posy,&posz,&refposx,&refposy,&refposz,&vx,&vy,&vz,&eigvects,&nnear,&inear,&rx,&ry,&rz,&modr,ebs,&lats,pbc,T,nu,ander,&TBparam);
	  ekin=3*(n-1)*kb*Tf/2;
	  if(i%nprint==0){
	    //H_MD and eigvects have now also been populated
	    erep=Erep(&modr);
	    etot=ebs+erep+ekin;
	    tmd=i*dt;
	    fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",tmd,Tf,ekin,ebs,erep,etot);
	    fprintf(mov,"%d\nC%d molecule\n",n,n);
	    for(j=0;j<n;j++){
	      fprintf(mov,"6  %f %f %f\n",posx.at(j),posy.at(j),posz.at(j));
	    }
	  } 
	  if(i*10%nmd==0){
	    std::cout << i*10/nmd << "0% completed" << std::endl;
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

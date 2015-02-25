#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "readinxyz.h"
#include "vectorfunctions.h"
#include "band_hamiltonian.h"
#include "functions.h"
#include "geometryinfo.h"

int main(int argc, char* argv[]){
	if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	
	// Turn verbose mode (hamiltonian routine) on/off
	bool v=0;
	// Read in types, 
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz);
	// Number of atoms
	int n=posx.size();
	std::cout << "n = " <<n <<std::endl;
	std::vector<double> posxnew(n), posynew(n), posznew(n);
	Eigen::MatrixXd modr(n,n);
	Eigen::MatrixXd rx(n,n);
	Eigen::MatrixXd ry(n,n);
	Eigen::MatrixXd rz(n,n);
	double ebs = 0,erep = 0,etot = 0;

	// Change if using PBCs
	double a =13.1;
	double b =10;
	double c =10;
	double rv = 2.98;

	Eigen::MatrixXd eigvects(4*n,4*n);

	// Calculate distances
	PbcGetAllDistances(&modr,&rx,&ry,&rz,&posx,&posy,&posz,a,b,c,rv);
	int kn = 5; // Number of kpts 
	std::vector<double> kvec(3);
	kvec.at(0) = 0;
	kvec.at(1) = 0;
	kvec.at(2) = 0;

	//k used for the band structure 
	//double k=0;
	// Starting TB	module: calculating energies
	std::vector<double> eigenvalaar (4*n);

	std::ofstream band;
	band.open ("band_structure6.dat");
	double K_MAX=M_PI/a;  //(a/n);
	double K_STEP=2*M_PI/(a*kn);        // 2*M_PI/((a/n)*50);
	double BZ = M_PI/(a/n);
	double kprim;
	//ebs=band_Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,&eigenvalaar,0,v);
	for(double k=0;k<=K_MAX;k=k+K_STEP){
	  kvec.at(0) = k;
	  std::cout << "k = " << k << std::endl;
	  ebs=ebs+(1.0/(double)kn)*band_Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,&eigenvalaar,&kvec,v);
		//H_MD and eigvects have now also been populated
		//band<<k<<"\t";
		for(int i=0;i<4*n;i++){

		  // BZ band 0
		  if (i/2-(double)i/2.0 == 0 && i < n) {             // If even band
			 kprim = fmod(i*K_MAX + k,BZ);
		  }
		  else if (i/2-(double)i/2.0 != 0 && i < n) {                                    // If odd band
			 kprim = fmod(i*K_MAX + (K_MAX-k),BZ);
		  }

		  // BZ band 1
		  else if (i/2-(double)i/2.0 == 0 && i < 2*n) {             // If even band
			 kprim = BZ - (i%n)*K_MAX - k;
		  }
		  else if (i/2-(double)i/2.0 != 0 && i < 2*n) {                                    // If odd band
			 kprim = BZ - (i%n)*K_MAX - (K_MAX - k);
		  }

		  // BZ band 2
		  if (i/2-(double)i/2.0 == 0 && i < 3*n) {             // If even band
			 kprim = fmod((i-2*n)*K_MAX + k,BZ);
		  }
		  else if (i/2-(double)i/2.0 != 0 && i < 3*n) {                                    // If odd band
			 kprim = fmod((i-2*n)*K_MAX + (K_MAX-k),BZ);
		  }

		  // Else
		  else if (i/2-(double)i/2.0 == 0 && i > n*2) {             // If even band
			 kprim = fmod(i*K_MAX + k,BZ);
		  }
		  else if (i/2-(double)i/2.0 != 0 && i > n*2) {                                    // If odd band
			 kprim = fmod(i*K_MAX + (K_MAX-k),BZ);
		  }

		  std::cout << "kprim = " << kprim << std::endl;
		  band<<kprim<<"\t";
		  band<<eigenvalaar[i]<<"\n";
		}
		//band<<"\n";
	}
	band.close();
	//H_MD and eigvects have now also been populated
	erep=Erep(&modr);
	// Determining erep works, however, we need to check if we're passing pointers or arrays to Erep()
	etot=ebs+erep;

	std::cout << "Ebs = " << ebs << std::endl;
	std::cout << "Erep = " << erep << std::endl;
	std::cout << "Etot = " << etot << std::endl;

	return 0;
}

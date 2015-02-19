#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "readinxyz.h"
#include "vectorfunctions.h"
#include "hamiltonian.h"
#include "functions.h"
#include "geometryinfo.h"
#include "ScaleGeom.h"

int main(int argc, char* argv[]){
	if (argc<3){std::cout<<"You should append three files to the main object!"<<std::endl;}
	if (argc!=4){std::cout<<"You should append two xyz files and one .dat file to the main!!"<<std::endl;}
	
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
	double ebs,erep,etot;

	double a =7.37854;
	double b =8.52;
	double c =13.416;
	double rv = 2.98;
	double bond = 1.42;

	std::string cellout = argv[2];
	std::string dataout = argv[3];
	std::ofstream cell(cellout.c_str());
	std::ofstream output(dataout.c_str());
	Eigen::MatrixXd eigvects(4*n,4*n);
	for (double scale=0.9;scale<=1.1;scale=scale+0.01){
	  std::cout << "Starting scale " << scale << std::endl;

	  // Scale cell size
	  ScaleGeom(n, scale, scale, 1, &posx, &posy, &posz, &posxnew, &posynew, &posznew);

	  // print .xyz file of cells
	  cell << n << "\n";
	  cell << "Scaled by " << scale << "\n";
	  for (int j=0;j<n;j++){
		 cell << "C\t" << posxnew.at(j) << "\t" << posynew.at(j) <<  "\t" << posznew.at(j) << "\n";
	  }
	  
	  // Calculate distances
	  PbcGetAllDistances(&modr,&rx,&ry,&rz,&posxnew,&posynew,&posznew,scale*a,scale*b,scale*c,rv);
	  
	  // Starting TB	module: calculating energies
	  ebs=Hamiltonian(n,&modr,&rx,&ry,&rz,&eigvects,v);
	  //H_MD and eigvects have now also been populated
	  erep=Erep(&modr);
	  // Determining erep works, however, we need to check if we're passing pointers or arrays to Erep()
	  etot=ebs+erep;
	  
	  output << scale*bond << "\t" << etot/(double)n << "\n";
	}
	
	return 0;
}

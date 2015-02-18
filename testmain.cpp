// Test for the distance calculations
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "geometryinfo.h"
#include "readinxyz.h"
#include "vectorfunctions.h"


int main(int argc, char* argv[]){
	if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
	// Read in types, 
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &posx, &posy, &posz);

	// Test distance function
	const int N=posx.size();
	Eigen::MatrixXd rx (N, N);
	Eigen::MatrixXd ry (N, N);
	Eigen::MatrixXd rz (N, N);
	Eigen::MatrixXd modr (N, N);
	std::vector<int>nnear(N);
	Eigen::MatrixXi inear = Eigen::MatrixXi::Constant(N, 10, -111);					// Random number to indicate that a value has been unassigned
	double rc=1.6, rv=1.6;
	
	double a=3.57, b=3.57, c=3.57;
	
	PbcGetAllDistances(&modr, &rx, &ry, &rz, &posx, &posy, &posz, a, b, c, rv);
	NearestNeighbours(&inear, &nnear, &modr, rv);
	
	
std::cout<<"rx: "<<std::endl;	
std::cout<<rx<<std::endl;
std::cout<<"ry: "<<std::endl;	
std::cout<<ry<<std::endl;
//std::cout<<"rz: "<<std::endl;	
//std::cout<<rz<<std::endl;
std::cout<<"modr: "<<std::endl;	
std::cout<<modr<<std::endl;
//std::cout<<"nnear: "<<std::endl;
//Print(&nnear);
//std::cout<<"inear: "<<std::endl;	
//std::cout<<inear<<std::endl;

return 0;
}

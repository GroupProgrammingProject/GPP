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
	std::vector<int> type;
	std::vector<double> posx, posy, posz;
	ReadInXYZ (argv[1], &type, &posx, &posy, &posz);

	// Test distance function
	const int N=type.size();
	Eigen::MatrixXd rx (N, N);
	Eigen::MatrixXd ry (N, N);
	Eigen::MatrixXd rz (N, N);
	Eigen::MatrixXd modr (N, N);
	std::vector<int>nnear(N);
	Eigen::MatrixXi inear = Eigen::MatrixXi::Constant(N, 5, -123456);					// Random number to indicate that a value has been unassigned
	GetAllDistances(&modr, &rx, &ry, &rz, &posx, &posy, &posz, 4);
	NearestNeighbours(&inear, &nnear, &modr, 4);
	
std::cout<<rx<<std::endl;
std::cout<<ry<<std::endl;
std::cout<<rz<<std::endl;
std::cout<<modr<<std::endl;
Print(&nnear);
std::cout<<inear<<std::endl;


//	std::vector<double> modr1(N*N), modr2(N*N), rx(N*N), ry(N*N), rz(N*N) ;
//	GetDistanceComponents(&rx,&ry, &rz, &posx, &posy, &posz);
//	GetAllDistances(&modr1, &rx, &ry, &rz, &posx, &posy, &posz);
//	GetAllDistances(&modr2, &rx, &ry, &rz, &posx, &posy, &posz, rc);

//	NearestNeighbours (&inear, &nnear, &modr1, 4.);
	
return 0;
}

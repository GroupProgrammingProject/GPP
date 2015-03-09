#include "../include/readinxyz.h"
// Function to initiate a map to translate element symbol to integer values corresponding to its atomic number
std::map<std::string, int> SetElementMap(){
	std::map<std::string, int> map;
	map["H"]=1;
	map["C"]=6;
	return map;
}

// Reads in 4 vectors from cell file: elements, x , y and z coords
void ReadInXYZ(char* filename, std::vector<double>* xvect, std::vector<double>* yvect, std::vector<double>* zvect, std::vector<double>* vxvect, std::vector<double>* vyvect, std::vector<double>* vzvect,std::vector<double>* lats, bool pbc, std::vector<bool>* velspec){
	std::ifstream infile(filename);
	// Initiate the map of elements (symbol-> atomic number)
	std::map<std::string, int> elmap=SetElementMap();
	// First line handled seperately and not to be stored
	std::string skipline1;
	int num, n;
	while (infile>>num){n=num;} //read in the number of atoms in the file
//	std::getline(infile, skipline1);
	// If using periodic boundary conditions
	if (pbc == 1) {
	  // Read in a, b and c lattice parameters
	  std::string name;
	  double a, b, c;
	  infile>>name>>a>>b>>c;
	  (*lats).at(0) = a;
	  (*lats).at(1) = b;
	  (*lats).at(2) = c;
	}
	else {
	  // Second line handled seperately and not to be stored
	  std::string skipline2;
	  std::getline(infile, skipline2);
	}
	// Store the molecule type and x, y, z positions and velocities
	double x, y, z, vx, vy, vz;
	std::string type;
	for(int i=0; i<n; i++)
	{
		std::string line;
		getline(infile, line);
		infile>>type>>x>>y>>z;
		xvect->push_back(x);
		yvect->push_back(y);
		zvect->push_back(z);
		if(line!=" "){
		vxvect->push_back(vx);
		vyvect->push_back(vy);
		vzvect->push_back(vz);
		}
	}
/*	while (infile>>type>>x>>y>>z){
		xvect->push_back(x);
		yvect->push_back(y);
		zvect->push_back(z);
	}*/
}

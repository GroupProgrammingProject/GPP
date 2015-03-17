#include "../include/readinxyz.h"
// Function to initiate a map to translate element symbol to integer values corresponding to its atomic number
std::map<std::string, int> SetElementMap(){
	std::map<std::string, int> map;
	map["H"]=1;
	map["C"]=6;
	return map;
}

// Reads in 4 vectors from cell file: elements; x , y and z coords; vx, vy and vz velocities; cell parameters; switch for PBCs; array of Booleans for input velocities
void ReadInXYZ(char* filename, std::vector<double>* xvect, std::vector<double>* yvect, std::vector<double>* zvect, std::vector<double>* vxvect, std::vector<double>* vyvect, std::vector<double>* vzvect,std::vector<double>* lats, bool pbc, std::vector<bool>* velspec){
	std::ifstream infile(filename);
	// Initiate the map of elements (symbol-> atomic number)
	std::map<std::string, int> elmap=SetElementMap();
	// First line handled seperately and not to be stored
	std::string skipline1, lstring;
	int n,i,j; //i,j dummy variables for 'for' loops
	std::getline(infile, skipline1);
	n=atoi(skipline1.c_str()); //read in number of atoms from first line
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
	std::string token; //store string input
	std::string delimiter=" "; //store delimiter
	size_t pos; //length of substring between delimiters
	for(i=0; i<n; i++){ //loop over number of lines specified at Line 1
		std::vector<std::string> input;
		std::string line;
		std::getline(infile, line);
		std::stringstream ss(line);
		while(ss >> token){input.push_back(token);}

/*		line=line+delimiter; //add a delimiter to the string so we can easily find end
		while(((pos=line.find(delimiter))<std::string::npos) && (line.length()>0)){
			token=line.substr(0,pos);
			if(token==" " || token==""){line.erase(0,pos+delimiter.length());} //get rid of token if we accidentally get a space or empty string 
			else{
				input.push_back (token);
				line.erase(0,pos+delimiter.length()); //delete part of string alreay read in
			}
		}
*//*		std::cout << "\nNumber of inputs on line " << i+2 << " =" << input.size() <<std::endl;
		for(j=0; j<input.size(); j++)
		{
			std::cout << input.at(j) << " ";//<< std::endl;
		}
*/	
		xvect->push_back(::atof(input.at(1).c_str())); //converts strings to doubles and adds the result to xvect vector
		yvect->push_back(::atof(input.at(2).c_str()));
		zvect->push_back(::atof(input.at(3).c_str()));
		if(input.size()>4) //if there are velocities on this line of the input file
		{
			vxvect->push_back(::atof(input.at(4).c_str()));
			vyvect->push_back(::atof(input.at(5).c_str()));
			vzvect->push_back(::atof(input.at(6).c_str()));
			velspec->push_back(1); //Boolean label tells us that there is an input velocity
		}
		else{ //if no input velocity is specified
			vxvect->push_back(0.0);
			vyvect->push_back(0.0);
			vzvect->push_back(0.0);
			velspec->push_back(0);
		}
	}
}

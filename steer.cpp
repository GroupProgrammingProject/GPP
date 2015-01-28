#include <iostream>
#include <fstream>
#include <vector>
#include "readinxyz.h"
#include "vectorfunctions.h"

int main(int argc, char* argv[]){
	std::vector<int> types;
	std::vector<double> X, Y, Z;
	std::vector<int>* typeptr = &types;
	std::vector<double>* Xptr = &X;
	std::vector<double>* Yptr = &Y;
	std::vector<double>* Zptr = &Z;
	
	//Read in values from specified .xyz file
	ReadInXYZ (argv[1], typeptr, Xptr, Yptr, Zptr);
	Print (typeptr);
	Print (Xptr);
	Print (Yptr);
	Print (Zptr);
	
	const char* filename="testz.txt";
	Print(Zptr, filename);
	
return 0;
}

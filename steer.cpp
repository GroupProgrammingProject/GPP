#include <iostream>
#include <fstream>
#include <vector>
#include "readinxyz.h"

// Function to print contents of an std::vector
template <typename T> void Print(std::vector<T>* vect){
	for (typename std::vector<T>::iterator it=vect->begin(); it!=vect->end(); ++it){
		std::cout<<*it<<"\t";
	}
	std::cout<<std::endl;
} 

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
return 0;
}

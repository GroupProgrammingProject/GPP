#ifndef VECTORFUNCTIONS_H
#define VECTORFUNCTIONS_H

// Some vector functions not included in the std::vector library
#include <iostream>
#include <fstream>
#include <vector>

template <typename T> void Print(std::vector<T>* vect);
template <typename T> void Print(std::vector<T>* vect, char* ofilename);

// Function to print contents of an std::vector to screen
template <typename T> void Print(std::vector<T>* vect){
	for (typename std::vector<T>::iterator it=vect->begin(); it!=vect->end(); ++it){
		std::cout<<*it<<"\t";
	}
	std::cout<<std::endl;
} 

// Function to print contents of an std::vector to file
template <typename T> void Print(std::vector<T>* vect, const char* filename){
	std::ofstream outfile(filename);
	for (typename std::vector<T>::iterator it=vect->begin(); it!=vect->end(); ++it){
		outfile<<*it<<"\t";
	}
} 

#endif

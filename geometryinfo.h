#ifndef GEOMETRYINFO_H
#define GEOMETRYINFO_H

#include <vector>
#include <cmath>
#include "vectorfunctions.h"
#include <Eigen/Dense>

// Function declarations

// Function to calculate one component of the distances between atoms
void DistanceComp(Eigen::MatrixXd* r, std::vector<double>* pos){
	int numatoms=pos->size();
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){(*r)(i, j)=pos->at(i)-pos->at(j);}
	}
}

// Function to calculate distance components between atoms
void GetDistanceComponents(Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	DistanceComp(rx, posx);
	DistanceComp(ry, posy);
	DistanceComp(rz, posz);
}


// Function to calculate the distances between atoms
void GetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	int numatoms=posx->size();
	GetDistanceComponents(rx, ry, rz, posx, posy, posz);
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){
		double distx=(*rx)(i, j);
		double disty=(*ry)(i, j);
		double distz=(*rz)(i, j);
		(*modr)(i, j)=sqrt(distx*distx+disty*disty+distz*distz);
		}
	}
}

// Function to calculate the distances between atoms and set all values above rc to zero
void GetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc){
	int numatoms=posx->size();
	GetDistanceComponents(rx, ry, rz, posx, posy, posz);
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){
		double distx=(*rx)(i, j);
		double disty=(*ry)(i, j);
		double distz=(*rz)(i, j);
		double modulus =sqrt(distx*distx+disty*disty+distz*distz);
		if (modulus>rc){(*modr)(i, j)=0;}
		else {(*modr)(i, j)= modulus;}
		}
	}
}


// Function that calculates the nearest neighbours of every atom
// Nearest neighbours are the neighbours within a radius rv
// inear is the vector in which the indices will be printed
void NearestNeighbours(Eigen::MatrixXi* inear, std::vector<int>* nnear, Eigen::MatrixXd* modr, double rv){
	int numatoms=nnear->size();
	for (int i=0; i<numatoms; i++){
		int nnearcounter=0;
		int k=0;
		for (int j =0; j<numatoms; j++){
			double dist = (*modr)(i, j);
			if (dist<rv && dist!=0){
				nnearcounter++;
				(*inear)(i, k)=j;
				k++;
			}
			else {}	
		}
		nnear->at(i)=nnearcounter;
	}
}


// Simple absolute value function
double Abs(double c){
	double result=c;
	if (c<0){result=-c;}
	else{;}
return result;
}

// Function to set all elements of 2 equal length vectors equal too each other
// equivalent to A=B for std::vector<double> objects
void SetEqual(std::vector<double>* a, std::vector<double>* b){
	int n=a->size();
	for (int i=0; i<n; i++){a->at(i)=b->at(i);}
}

// Function that decides whether the Nearest neighbours lists need to be recalculated between MD steps
bool RecalculateNearestNeighbours(std::vector<double>* refposx, std::vector<double>* refposy, std::vector<double>* refposz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc, double rv){
	int numatoms=posx->size();
	bool recalc=0;
	std::vector<double> disp(numatoms);
	for (int i=0; i<numatoms; i++){
		double xdisp=posx->at(i)-refposx->at(i);
		double ydisp=posy->at(i)-refposy->at(i);
		double zdisp=posz->at(i)-refposz->at(i);
		disp.at(i)=sqrt(xdisp*xdisp+ydisp*ydisp+zdisp*zdisp);
	}
	double dmax=*(max_element(disp.begin(), disp.end()));
std::cout<<dmax<<std::endl;
	if (dmax>0.4*Abs(rc-rv)){
		SetEqual(posx, refposx);
		SetEqual(posy, refposy);
		SetEqual(posz, refposz);
		recalc=1;
	}
	else{;}
return recalc;
}


#endif

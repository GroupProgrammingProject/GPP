#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "Gethijab.h"
#include "functions.h"

// Takes 1 int and 4 vector arguments: n, type, posx, posy, posz and returns band structure energy Ebs
double Hamiltonian(int n, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* eigvects){
  
  int i, j, a, b;                                         // i,j loop over atoms; a,b loop over orbitals
  double r,rx,ry,rz;                                 // d is a 3d array of atom pair's connecting vector (varys within ij loop)
  std::vector<double> d(3);
  double sr,hijab;                                        // hijab: unscaled matrix element, sr: scaling function value at r
  double Ebs=0;														 // Band structure energy to be returned
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  
  MatrixXd Hijab(4*n,4*n);											 // Hamiltonian matrix for diagonalization

  for (i=0;i<n;i++) {                                     // Cycle through atoms i
    for (j=0;j<=i;j++) {                                  // Cycle through atoms j (note upper triangular matrix)
      rx    = (*posx).at(i)-(*posx).at(j);
      ry    = (*posy).at(i)-(*posy).at(j);                // Cartesian elements of vector r[i] - r[j]
      rz    = (*posz).at(i)-(*posz).at(j);
      r = sqrt(pow(rx,2)+ pow(ry,2) + pow(rz,2));	       // Length |r[i] - r[j]|
 		rx=rx/r;
		ry=ry/r;
		rz=rz/r;
      d.at(0)=rx; d.at(1)=ry; d.at(2)=rz;
		  	if (i == j) {sr = 1;}                               // Don't apply scaling function if i=j
      	else {sr    = s(r);}                  // Scaling parameter
      	for (a=0;a<4;a++) {                                 // Cycle through orbitals of atom i
				for (b=0;b<4;b++) {                               // Cycle through orbitals of atom i
	  				if (sr == 0) {hijab = 0;}                       // If scaling function gives 0, no need to calc hijab
	  				else {hijab = Gethijab(i,j,a,b,&d);} // Hamiltonian elements of ij interaction
	  				Hijab(4*i+a,4*j+b)     = sr*hijab;              // Scale hijab and populate matrix Hijab
	  				Hijab(4*j+b,4*i+a)     = sr*hijab;              // Scale hijab and populate matrix Hijab
				}                                                 // End loop over b
      	}                                                   // End loop over a
    	}                                                     // End loop over j
  	}                                                       // End loop over i


  Eigen::SelfAdjointEigenSolver<MatrixXd> es(Hijab);         // Compute eigenvectors and eigenvalues

  typedef std::pair<double,int> pair;
  std::vector<pair> eigvalarr;
  for(i=0;i<4*n;i++){eigvalarr.push_back(pair(es.eigenvalues()[i],i));}		//reads in eigenvalues and their index into eigvalarr
  std::sort(eigvalarr.begin(),eigvalarr.end());										//sorts eigenvalues
  for (i=0;i<2*n;i++) {Ebs = Ebs + 2*eigvalarr.at(i).first;}           		// Fill lowest eigenstates with 2 electrons and sum energies of filled states

  // Filling up eigvects, doubling the vectors of the states that are occupied
  for(i=0;i<2*n;i++){
  		int ind=eigvalarr.at(i).second;
	  for(j=0;j<4*n;j++){
			(*eigvects).at(2*i*4*n+j)=es.eigenvectors().row(j).col(ind).value();			//reads in eigenvectors for occupied states only
			(*eigvects).at((2*i+1)*4*n+j)=es.eigenvectors().row(j).col(ind).value();	//twice, as each state is "doubly-occupied"
	  }
  }

// Uncomment any section to print output for testing

/* 	std::cout << Hijab << std::endl;		//print out Hamiltonian
	  std::cout << "After sorting" << std::endl;
  for(i=0;i<4*n;i++){std::cout << "eigvalarr no. " << eigvalarr.at(i).second << " is " << eigvalarr.at(i).first << std::endl;}
	
    std::cout << "Eigenvector matrix" << std::endl;
  std::cout << es.eigenvectors() << std::endl;									//untouched eigenvector matrix

  std::cout << "Eigenvectors after sorting and filling only occupied states:" << std::endl;
  for(i=0;i<4*n;i++){
		for(j=0;j<4*n;j++){
			std::cout << (*eigvects).at(i*4*n+j) << "\t\t";
		}
	std::cout << std::endl;
	}
*/

  return Ebs;
}

#endif

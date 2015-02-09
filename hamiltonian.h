#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "Gethijab.h"
#include "functions.h"

// Need to make number of orbitals (currently m=4) a variable since needs to account for total number of electrons
// WE DON'T Consider systems with odd numbers of electrons

// Takes 1 int and 4 vector arguments: n, type, posx, posy, posz and returns band structure energy Ebs
double Hamiltonian(int n, std::vector<int>* type, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* eigvects){
  
  int i, j, a, b;                                         // i,j loop over atoms; a,b loop over orbitals
  double d[3],r,rx,ry,rz;                                 // d is a 3d array of atom pair's connecting vector (varys within ij loop)
  double sr,hijab;                                        // hijab: unscaled matrix element, sr: scaling function value at r
  double Ebs=0;														 // Band structure energy to be returned
  int  typei, typej;                                      // Atomic numbers of atoms i and j
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> VectorXd;
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixComplex;
  
  MatrixXd Hijab(4*n,4*n);											 // Hamiltonian matrix for diagonalization

  for (i=0;i<n;i++) {                                     // Cycle through atoms i
    for (j=0;j<=i;j++) {                                  // Cycle through atoms j (note upper triangular matrix)
      typei = (*type).at(i);                              // Extract element of atom i
      typej = (*type).at(j);      
      rx    = (*posx).at(i)-(*posx).at(j);
      ry    = (*posy).at(i)-(*posy).at(j);                // Cartesian elements of vector r[i] - r[j]
      rz    = (*posz).at(i)-(*posz).at(j);
      d[0]=rx; d[1]=ry; d[2]=rz;
      r = sqrt(pow(rx,2)+ pow(ry,2) + pow(rz,2));	       // Length |r[i] - r[j]|
      	if (r == 0) {sr = 1;}                               // Don't apply scaling function if i=j
      	else {sr    = s(r, typei, typej);}                  // Scaling parameter
      	for (a=0;a<4;a++) {                                 // Cycle through orbitals of atom i
				for (b=0;b<4;b++) {                               // Cycle through orbitals of atom i
	  				if (sr == 0) {hijab = 0;}                       // If scaling function gives 0, no need to calc hijab
	  				else {hijab = Gethijab(i,j,a,b,d,typei,typej);} // Hamiltonian elements of ij interaction
	  				Hijab(4*i+a,4*j+b)     = sr*hijab;              // Scale hijab and populate matrix Hijab
	  				Hijab(4*j+b,4*i+a)     = sr*hijab;              // Scale hijab and populate matrix Hijab
				}                                                 // End loop over b
      	}                                                   // End loop over a
    	}                                                     // End loop over j
  	}                                                       // End loop over i

 	std::cout << Hijab << std::endl;

  Eigen::SelfAdjointEigenSolver<MatrixXd> es(Hijab);         // Compute eigenvectors and eigenvalues
//	Eigen::ComplexEigenSolver<MatrixXd> ces(Hijab);
//  VectorXd eigvals = ces.eigenvalues();                     // Retrieve Eigenvalues

  typedef std::pair<double,int> pair;
  std::vector<pair> eigvalarr;
  for(i=0;i<4*n;i++){
  		eigvalarr.push_back(pair(es.eigenvalues()[i],i));
//		eigvalarr.push_back(pair(real(ces.eigenvalues()[i]),i));
//  		std::cout << "eigvalarr no. " << eigvalarr.at(i).second << " is " << eigvalarr.at(i).first << std::endl;
  }

  std::sort(eigvalarr.begin(),eigvalarr.end());
//  std::cout << "After sorting" << std::endl;
//  for(i=0;i<4*n;i++){std::cout << "eigvalarr no. " << eigvalarr.at(i).second << " is " << eigvalarr.at(i).first << std::endl;}

//	std::cout << "Difference of degenerate values is" << eigvalarr.at(3).first-eigvalarr.at(4).first << std::endl;

  for (i=0;i<n;i++) {Ebs = Ebs + 2*eigvalarr.at(i).first;}           // Fill lowest eigenstates with 2 electrons and sum energies of filled states

//  MatrixComplex eigvectors = ces.eigenvectors();					// Retrieve Eigenvectors

//  std::cout << "Eigenvector matrix" << std::endl;
//  std::cout << es.eigenvectors() << std::endl;
//	std::cout << ces.eigenvectors() << std::endl;
  for(i=0;i<4*n;i++){
  		int ind=eigvalarr.at(i).second;
//		std::cout << "index of current eigval is " << ind << std::endl;
	  for(j=0;j<4*n;j++){
			(*eigvects).at(i*4*n+j)=es.eigenvectors().row(j).col(ind).value();
//real(eigvectors(j,ind));
//			std::cout << (*eigvects).at(j*4*n+i) << std::endl;
	  }
  }

/*
  std::cout << "Sorted eigenvectors in matrix:" << std::endl;
  for(i=0;i<4*n;i++){
		for(j=0;j<4*n;j++){
			std::cout << (*eigvects).at(j*4*n+i) << "		";
		}
	std::cout << std::endl;
	}
*/


  return Ebs;
}

#endif

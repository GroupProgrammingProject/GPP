#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "Gethijab.h"
#include "S.h"

// Takes 1 int and 4 vector arguments: n, type, posx, posy, posz
void Hamiltonian(int n, std::vector<int>* type, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
  int i, j, a, b;                                         // i,j loop over atoms; a,b loop over orbitals
  double d[3], r;                                         // 3d array of atom pair's connecting vector (varys within ij loop)
  double sr;                                              // Value of scaling functn for each ij pair
  double hijab;                                           // Matrix element to be rewritten during loop
  int  typei, typej;                                      // Atomic numbers of atoms i and j
  std::vector<double> H_MD(pow(4*n,2));                   // A matrix to be populated with interaction parameters for MD
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  MatrixXd Hijab(4*n,4*n);

  for (i=0;i<n;i++) {                                     // Cycle through atoms i
    for (j=0;j<=i;j++) {                                  // Cycle through atoms j (note upper triangular matrix)
      typei = (*type).at(i);                              // Extract element of atom i
      typej = (*type).at(j);      
      d[0]  = (*posx).at(i)-(*posx).at(j);
      d[1]  = (*posy).at(i)-(*posy).at(j);                // Cartesian elements of vector r[i] - r[j]
      d[2]  = (*posz).at(i)-(*posz).at(j);
      r     = sqrt(pow(rx,2)+ pow(ry,2) + pow(rz,2));     // Length |r[i] - r[j]|
      sr    = S(r, typei, typej);                         // Scaling parameter
      for (a=0;a<4;a++) {                                 // Cycle through orbitals of atom i
	for (b=0;b<4;b++) {                               // Cycle through orbitals of atom i
	  if (sr == 0) {hijab = 0;}                       // If scaling function gives 0, no need to calc hijab
	  else {hijab = Gethijab(i,j,a,b,d,typei,typej);} // Hamiltonian elements of ij interaction
	  H_MD.at((4*i+a)*4*j+b) = hijab;                 // Vector of interactions to pass to MD
	  H_MD.at((4*j+a)*4*i+b) = hijab;                 // Vector of interactions to pass to MD
	  Hijab(4*i+a,4*j+b)     = sr*hijab;              // Scale hijab and populate matrix Hijab
	  Hijab(4*j+a,4*i+b)     = sr*hijab;              // Scale hijab and populate matrix Hijab
	}                                                 // End loop over b
      }                                                   // End loop over a
    }                                                     // End loop over j
  }                                                       // End loop over i

  std::cout << Hijab << std::endl;                        // Print Hamiltonian Matrix

}

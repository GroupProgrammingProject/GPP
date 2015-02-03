#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "Gethijab.h"
#include "functions2.h"

// Need to make number of orbitals (currently m=4) a variable since needs to account for total number of electrons
// WE DON'T Consider systems with odd numbers of electrons

// Takes 1 int and 4 vector arguments: n, type, posx, posy, posz
double Hamiltonian(int n, std::vector<int>* type, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
  int i, j, a, b;                                         // i,j loop over atoms; a,b loop over orbitals
  double d[3], r;                                         // 3d array of atom pair's connecting vector (varys within ij loop)
  double rx, ry, rz;
  double sr;                                              // Value of scaling functn for each ij pair
  double hijab;                                           // Matrix element to be rewritten during loop
  int  typei, typej;                                      // Atomic numbers of atoms i and j
  std::vector<double> H_MD(pow(4*n,2));                   // A matrix to be populated with interaction parameters for MD
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  MatrixXd Hijab(4*n,4*n);
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> VectorXd;

  for (i=0;i<n;i++) {                                     // Cycle through atoms i
    for (j=0;j<=i;j++) {                                  // Cycle through atoms j (note upper triangular matrix)
      typei = (*type).at(i);                              // Extract element of atom i
      typej = (*type).at(j);      
      rx    = (*posx).at(i)-(*posx).at(j);
      ry    = (*posy).at(i)-(*posy).at(j);                // Cartesian elements of vector r[i] - r[j]
      rz    = (*posz).at(i)-(*posz).at(j);
      d[0]  = rx;
      d[1]  = ry;
      d[2]  = rz;
      r     = sqrt(pow(rx,2)+ pow(ry,2) + pow(rz,2));     // Length |r[i] - r[j]|
      if (r == 0) {sr = 1;}                               // Don't apply scaling function if i=j
      else {sr    = s(r, typei, typej);}                  // Scaling parameter
      for (a=0;a<4;a++) {                                 // Cycle through orbitals of atom i
	for (b=0;b<4;b++) {                               // Cycle through orbitals of atom i
	  if (sr == 0) {hijab = 0;}                       // If scaling function gives 0, no need to calc hijab
	  else {hijab = Gethijab(i,j,a,b,d,typei,typej);} // Hamiltonian elements of ij interaction
	  H_MD.at((4*i+a)*4*j+b) = hijab;                 // Vector of interactions to pass to MD
	  H_MD.at((4*j+b)*4*i+a) = hijab;                 // Vector of interactions to pass to MD
	  Hijab(4*i+a,4*j+b)     = sr*hijab;              // Scale hijab and populate matrix Hijab
	  Hijab(4*j+b,4*i+a)     = sr*hijab;              // Scale hijab and populate matrix Hijab
	}                                                 // End loop over b
      }                                                   // End loop over a
    }                                                     // End loop over j
  }                                                       // End loop over i

  std::cout << Hijab << std::endl;                        // Print Hamiltonian Matrix
  VectorXd eigvals = Hijab.eigenvalues();                 // Evaluate Eigenvalues
  double eigvalarr[4*n];                                  // Create array of Eigenvalues
  for (i=0;i<4*n;i++) {
    eigvalarr[i] = real(eigvals(i,0));}
  std::sort(eigvalarr,eigvalarr+4*n);                     // Sort array of eigenvalues
  double Ebs=0;                                           // Introduce band structure energy
  for (i=0;i<n;i++) {                                     // Fill each eigenstate with 2 electrons
    Ebs = Ebs + 2*eigvalarr[i];}                          // and sum energies of filled states

  std::cout << "\nEigenvalues of Hijab:\n" << std::endl;  // Output Eigenvalues
  for (i=0;i<4*n;i++) {
    std::cout << eigvalarr[i] << std::endl;}
  std::cout << "\nEbs = " << Ebs << "eV" << std::endl;    // Output Band Structure Energy

  // NOTE may be computing diagonalisation for a second time at this point
  // This is clearly unnecessary
  Eigen::ComplexEigenSolver<MatrixXd> ces(Hijab);         // Compute eigenvectors

  // Uncomment to view eigenvectors
  /*std::cout << "\nEigenvectors of Hijab are:\n" << std::endl;
  
  for (i=0;i<4*n;i++) {
    std::cout << "\nFor eigenvalue " << eigvals(i,0) << ":\n" << std::endl;    
    std::cout << ces.eigenvectors().col(i) << std::endl;
    }*/

  return Ebs;

}

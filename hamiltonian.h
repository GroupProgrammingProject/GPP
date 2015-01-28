#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "Gethijab.h"
#include "S.h"

/*double Gethijab(int i, int j, int a, int b, double* d, double r, int typei, int typej);
  double S(double r);*/

int Test(int n, std::vector<double>* posx) {
  return 1;
}

// Takes 1 int and 4 vector arguments: n, type, posx, posy, posz
void Hamiltonian(int n, std::vector<int>* type, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
  int i, j, a, b;                                       // i,j loop over atoms; a,b loop over orbitals
  double d[3], r;                                       // 3d array of atom pair's connecting vector (varys within ij loop)
  double rx, ry, rz, sr;
  double hijab;
  int  typei, typej;
  std::vector<double> H_MD(pow(4*n,2));                 // A matrix to be populated with interaction parameters for MD
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  MatrixXd Hijab(4*n,4*n);

  std::cout << "Checkpoint A" << std::endl;

  for (i=0;i<n;i++) {                                   // Cycle through atoms i
    std::cout << "i = " << i << std::endl;
    for (j=0;j<i;j++) {                                 // Cycle through atoms j (note upper triangular matrix)
      std::cout << "j = " << j << std::endl;
      typei = (*type).at(i);
      typej = (*type).at(j);      
      d[0]  = (*posx).at(i)-(*posx).at(j);
      d[1]  = (*posy).at(i)-(*posy).at(j);              // Cartesian elements of vector r[i] - r[j]
      d[2]  = (*posz).at(i)-(*posz).at(j);
      r     = sqrt(pow(rx,2)+ pow(ry,2) + pow(rz,2));   // Length |r[i] - r[j]|
      sr    = S(r);                                     // Scaling parameter
      for (a=0;a<4;a++) {
	for (b=0;b<4;b++) {
	  std::cout << "hijab = h" << i << j << a << b << std::endl;
	  hijab = Gethijab(i,j,a,b,d,r,typei,typej);    // Hamiltonian elements of ij interaction
	  H_MD.at((4*i+a)*4*j+b) = hijab;               // Vector of interactions to pass to MD
	  H_MD.at((4*j+a)*4*i+b) = hijab;               // Vector of interactions to pass to MD
	  std::cout << "hijab = " << hijab << std::endl;
	  Hijab(4*i+a,4*j+b)     = sr*hijab;            // Scale hijab and populate matrix Hijab
	  Hijab(4*j+a,4*i+b)     = sr*hijab;            // Scale hijab and populate matrix Hijab
	}                                               // End loop over b
      }                                                 // End loop over a
    }                                                   // End loop over i
  }                                                     // End loop over j

  std::cout << "Checkpoint B" << std::endl;

  std::cout << Hijab << std::endl;

}

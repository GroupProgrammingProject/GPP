#include <iostream>
#include <vector>
#include <Eigen/Dense>

double get_hijab(int a, int b, double* d);
int tot_orbitals(std::string atom);
double s(double r);

// Takes 1 int and 4 vector arguments: n, type, posx, posy, posz
void hamiltonian(int n, std::string* type, double* posx, double* posy, double* posz){
  int i, j, a, b;                          // i,j loop over atoms; a,b loop over orbitals
  double d[3];                             // 3d array of atom pair's connecting vector (varys within ij loop)
  std::vector<double> H_MD;                // A matrix to be populated with interaction parameters for MD
  typedef Eigen::Matrix<std::complex<double>, Dynamic, Dynamic > Matrix;
  Matrix Hijab;

  int row = 0;                             // Keep track of row in Hijab
  for (i=0:i<n;i++) {                      // Cycle through atoms i
    int col=0, col0=0;                     // Keep track of column in Hijab, col0 is origin of mini-matrix
    for (j=0:j<i;j++) {                    // Cycle through atoms j (note upper triangular matrix)
      int mi = 4; //int mi = tot_orbitals(type[i]); // Number of orbitals of atom i
      int mj = 4; //int mj = tot_orbitals(type[j]); // Number of orbitals of atom j
      for (a=0:a<mi;a++) {row++;
	col = col0;                         // returned to begining column of mini-matrix, update col
	for (b=0:b<mj;b++) {col++;
	  d = [posx[i]-posx[j],posy[i]-posy[j],posz[i]-posz[j]];  // Vector r[i] - r[j]
	  H_MD[] = get_hijab(a,b,d);        // Hamiltonian elements of ij interaction to pass to MD
	  Hijab(row,col);
	} col0 = col;                       // End loop over b, update zero column
      } row0 = row;                         // End loop over a, update zero row
    }} // End loop over i,j
}

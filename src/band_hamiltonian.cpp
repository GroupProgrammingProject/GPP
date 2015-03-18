#include "../include/band_hamiltonian.h"
// Takes scalar product between two vectors a and b
double dot(std::vector<double>* a, std::vector<double>* b) {
  if ((*a).size() != (*b).size()) {
	 std::cout << "Error: Cannot take scalar product of vectors of unequal length!" << std::endl;
  }
  else {
	 double product = 0;
	 for (int i=0;i<(*a).size();i++) {
		product = product + (*a).at(i)*(*b).at(i);
	 }
	 return product;
  }
}

// Takes 1 int and 4 vector arguments: n, type, posx, posy, posz and returns band structure energy Ebs
double band_Hamiltonian(int n,int norbs,std::vector<double>* TBparam,Eigen::MatrixXd* modr, Eigen::MatrixXd* rx,Eigen::MatrixXd* ry,Eigen::MatrixXd* rz,Eigen::MatrixXcd* eigvects,std::vector<double>* eigenvalaar, std::vector<double>* k, bool verbose){
  
  int i, j, a, b;                                          // i,j loop over atoms; a,b loop over orbitals
  double r,rxloc, ryloc, rzloc;                            // temporary parameters to store interatomic distances
  std::vector<double> d(3);                                // d is a 3d array of atom pair's connecting vector (varys within ij loop)
  double sr,hijab;                                         // hijab: unscaled matrix element, sr: scaling function value at r
  double Ebs=0;						                          // Band structure energy to be returned
  Eigen::MatrixXcd Hijab(norbs*n,norbs*n);			        // Hamiltonian matrix for diagonalization (complex this time)
  double kr;                                               // scalar product of k and r vectors (to avoid recalculating)

  for (i=0;i<n;i++) {                                      // Cycle through atoms i
    for (j=0;j<=i;j++) {                                   // Cycle through atoms j (note upper triangular matrix)
		r = (*modr)(i, j);
		if (i == j) {sr = 1; kr = 0;}                        // Don't apply scaling function if i=j
		else if (r < 1e-5) {sr = 0; kr = 0;}                 // If atoms have no interaction distance
		else {                                               // If atoms have finite interaction distance
		  sr    = s(r);                                      // Scaling parameter
		  rxloc = (*rx)(i, j)/r;
		  ryloc = (*ry)(i, j)/r;
		  rzloc = (*rz)(i, j)/r;
		  d.at(0)=rxloc; d.at(1)=ryloc; d.at(2)=rzloc;
		  kr = dot(k,&d)*r;
		}
      	for (a=0;a<norbs;a++) {                               // Cycle through orbitals of atom 
				for (b=0;b<norbs;b++) {                            // Cycle through orbitals of atom i
	  				if (sr == 0) {hijab = 0;}                       // If scaling function gives 0, no need to calc hijab
	  				else {hijab = Gethijab(i,j,a,b,&d,TBparam);}    // Hamiltonian elements of ij interaction
	  				Hijab(norbs*i+a,norbs*j+b) = exp(std::complex<double>(0,1)*kr)*sr*hijab;       // Scale hijab and populate matrix Hijab
	  				Hijab(norbs*j+b,norbs*i+a) = exp(-std::complex<double>(0,1)*kr)*sr*hijab;      // Scale hijab and populate matrix Hijab
				}                                              // End loop over b
      	}                                                 // End loop over a
    	}                                                    // End loop over j
  	}                                                       // End loop over i

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hijab);       // Compute eigenvectors and eigenvalues

  typedef std::pair<double,int> pair;                              // Create pair class for storing matched indicies

  std::vector<pair> eigvalarr;
  for(i=0;i<norbs*n;i++){eigvalarr.push_back(pair(es.eigenvalues()[i],i));}		// reads in eigenvalues and their index into eigvalarr
  std::sort(eigvalarr.begin(),eigvalarr.end());										      // sorts eigenvalues
  for (i=0;i<2*n;i++) {Ebs = Ebs + 2*eigvalarr.at(i).first;}       	            // Fill lowest eigenstates with 2 electrons and sum energies of filled states

  for(i=0;i<norbs*n;i++){(*eigenvalaar).at(i) = eigvalarr.at(i).first;}

  for(i=0;i<2*n;i++){
  		int ind=eigvalarr.at(i).second;
	  for(j=0;j<norbs*n;j++){
			(*eigvects)(2*i,j)=es.eigenvectors().row(j).col(ind).value();			//reads in eigenvectors for occupied states only
			(*eigvects)(2*i+1,j)=es.eigenvectors().row(j).col(ind).value();	   //twice, as each state is "doubly-occupied"
	  }
  }

  // Output to terminal, set verbose == 1 to print
  if (verbose == 1) {
	 std::cout << "\nReal Hijab:"  << std::endl;                                           //print out Hamiltonian
	 std::cout << Hijab.real() << std::endl;		                                          //print out Hamiltonian
	 std::cout << "\nImaginary Hijab:"  << std::endl;                                      //print out Hamiltonian
	 std::cout << Hijab.imag() << std::endl;		                                          //print out Hamiltonian
	 std::cout << "After sorting" << std::endl;
	 for(i=0;i<norbs*n;i++){std::cout << "eigvalarr no. " << eigvalarr.at(i).second << " is " << eigvalarr.at(i).first << std::endl;}
	 
	 std::cout << "Eigenvector matrix" << std::endl;
	 std::cout << es.eigenvectors() << std::endl;									               //untouched eigenvector matrix
	 
	 std::cout << "Eigenvectors after sorting and filling only occupied states:" << std::endl;
	 std::cout << "\nReal eigenvectors:"  << std::endl;                                    //print out Hamiltonian
	 std::cout << (*eigvects).real() << std::endl;
	 std::cout << "\nImaginary eigenvectors:"  << std::endl;                               //print out Hamiltonian
	 std::cout << (*eigvects).imag() << std::endl;
  } // End of verbose mode

  return Ebs;
}



#include "../include/hamiltonian.h"

double Hamiltonian(int n,std::vector<double>* TBparam, Eigen::MatrixXd* modr, Eigen::MatrixXd* rx,Eigen::MatrixXd* ry,Eigen::MatrixXd* rz,Eigen::MatrixXd* eigvects, bool verbose){
  int i, j, a, b;                                          // i,j loop over atoms; a,b loop over orbitals
  double r,rxloc, ryloc, rzloc;                            // temporary parameters to store interatomic distances
  std::vector<double> d(3);                                // d is a 3d array of atom pair's connecting vector (varys within ij loop)
  double sr,hijab;                                         // hijab: unscaled matrix element, sr: scaling function value at r
  double Ebs=0;										   			  // Band structure energy to be returned
  Eigen::MatrixXd Hijab(4*n,4*n);			   							  // Hamiltonian matrix for diagonalization

  for (i=0;i<n;i++) {                                      // Cycle through atoms i
    for (j=0;j<=i;j++) {                                   // Cycle through atoms j (note upper triangular matrix)
		r = (*modr)(i, j);
		if (i == j) {sr = 1;}                                // Don't apply scaling function if i=j
		else if (r < 1e-5) {sr = 0;}                         // If atoms have no interaction distance
		else {                                               // If atoms have finite interaction distance
		  sr    = s(r);                                      // Scaling parameter
		  rxloc = (*rx)(i, j)/r;
		  ryloc = (*ry)(i, j)/r;
		  rzloc = (*rz)(i, j)/r;
		  d.at(0)=rxloc; d.at(1)=ryloc; d.at(2)=rzloc;
		}
		for (a=0;a<4;a++) {                               // Cycle through orbitals of atom i
		  for (b=0;b<4;b++) {                            // Cycle through orbitals of atom i
			 if (sr == 0) {hijab = 0;}                   // If scaling function gives 0, no need to calc hijab
			 else {hijab = Gethijab(i,j,a,b,&d,TBparam);}        // Hamiltonian elements of ij interaction
			 Hijab(4*i+a,4*j+b)     = sr*hijab;          // Scale hijab and populate matrix Hijab
			 //	  				Hijab(4*j+b,4*i+a)     = sr*hijab;          // Scale hijab and populate matrix Hijab
		  }                                              // End loop over b
		}                                                 // End loop over a
	 }                                                    // End loop over j
  }                                                       // End loop over i

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hijab);       // Compute eigenvectors and eigenvalues

  typedef std::pair<double,int> pair;                                      // Create pair class for storing matched indicies
  std::vector<pair> eigvalarr;                                             
  for(i=0;i<4*n;i++){eigvalarr.push_back(pair(es.eigenvalues()[i],i));}		// reads in eigenvalues and their index into eigvalarr
  std::sort(eigvalarr.begin(),eigvalarr.end());										// sorts eigenvalues
  for (i=0;i<2*n;i++) {Ebs = Ebs + 2*eigvalarr.at(i).first;}           		// Fill lowest eigenstates with 2 electrons and sum energies of filled states

  for(i=0;i<2*n;i++){
  		int ind=eigvalarr.at(i).second;
	  for(j=0;j<4*n;j++){
			(*eigvects)(2*i,j)=es.eigenvectors().row(j).col(ind).value();			//reads in eigenvectors for occupied states only
			(*eigvects)(2*i+1,j)=es.eigenvectors().row(j).col(ind).value();	//twice, as each state is "doubly-occupied"
	  }
  }

  // Output to terminal, set verbose == 1 to print
  if (verbose == 1) {
	 std::cout << Hijab << std::endl;		                                          //print out Hamiltonian
	 std::cout << "After sorting" << std::endl;
	 for(i=0;i<4*n;i++){std::cout << "eigvalarr no. " << eigvalarr.at(i).second << " is " << eigvalarr.at(i).first << std::endl;}
	 
	 std::cout << "Eigenvector matrix" << std::endl;
	 std::cout << es.eigenvectors() << std::endl;									         //untouched eigenvector matrix
	 
	 std::cout << "Eigenvectors after sorting and filling only occupied states:" << std::endl;
	 for(i=0;i<4*n;i++){
		for(j=0;j<4*n;j++){
		  std::cout << (*eigvects)(j,i) << "\t\t";
		}
		std::cout << std::endl;
	 }
  } // End of verbose mode


  return Ebs;
}



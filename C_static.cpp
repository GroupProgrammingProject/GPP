#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "hamiltonian.h"
#include "functions.h"

int main() {
  int n=3;
  std::vector<double> posx(n);
  std::vector<double> posy(n);
  std::vector<double> posz(n);
  std::vector<double> eigvects(16*n*n);
  for(int i=0;i<n;i++){
	posy.at(i)=0.0;
	posz.at(i)=0.0;
}

  std::ofstream output("C_bondenergy.dat");
  for (double bond=1;bond<=2;bond=bond+0.01){
	for(int j=0;j<n;j++){
		posx.at(j)=j*bond;
	}
    double ebs = Hamiltonian(n, &posx, &posy, &posz, &eigvects);
    double erep = Erep(&posx, &posy, &posz);  
    double etot = ebs + erep;
    output << bond << "\t" << etot/n << "\n";
  }

  return 0;
}

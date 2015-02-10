#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "hamiltonian.h"
#include "functions.h"

int main() {
  int n=8;
  std::vector<double> posx(n);
  std::vector<double> posy(n);
  std::vector<double> posz(n);
  std::vector<double> eigvects(16*n*n);
  posx.at(0) = 0.0000;
  posy.at(0) = 0.0000;
  posy.at(1) = 0.0000;
  posz.at(0) = 0.0000;
  posz.at(1) = 0.0000;
  posy.at(2) = 0.0000;
  posz.at(2) = 0.0000;
  posy.at(3) = 0.0000;
  posz.at(3) = 0.0000;
  posy.at(4) = 0.0000;
  posz.at(4) = 0.0000;
  posy.at(5) = 0.0000;
  posz.at(5) = 0.0000;
  posy.at(6) = 0.0000;
  posz.at(6) = 0.0000;
  posy.at(7) = 0.0000;
  posz.at(7) = 0.0000;

  std::ofstream output("C8bondenergy.dat");
  for (double bond=1;bond<=2;bond=bond+0.01){
    posx.at(1) = bond;
    posx.at(2) = 2*bond;
    posx.at(3) = 3*bond;
    posx.at(4) = 4*bond;
    posx.at(5) = 5*bond;
    posx.at(6) = 6*bond;
    posx.at(7) = 7*bond;

    double ebs = Hamiltonian(n, &posx, &posy, &posz, &eigvects);
    double erep = Erep(&posx, &posy, &posz);  
    double etot = ebs + erep;
    output << bond << "\t" << etot << "\n";
  }

  return 0;
}

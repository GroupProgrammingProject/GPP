#include <iostream>
#include <vector>
#include <cmath>
#include "hamiltonian.h"
#include "functions.h"

int main() {
  std::vector<int> type(2);
  std::vector<double> posx(2);
  std::vector<double> posy(2);
  std::vector<double> posz(2);
  int i,j,n = 2;
  std::vector<double> H_MD(16*n*n);
  std::vector<double> eigvects(16*n*n);

  type.at(0) = 6;
  type.at(1) = 6;
  posx.at(0) = 0.0000;
  posx.at(1) = 1.2400;
  posy.at(0) = 0.0000;
  posy.at(1) = 0.0000;
  posz.at(0) = 0.0000;
  posz.at(1) = 0.0000;

  double ebs = Hamiltonian(n, &type, &posx, &posy, &posz, &H_MD, &eigvects);
  double erep = Erep(type, posx, posy, posz);  
  double etot = ebs + erep;

  std::cout << "Ebs = " << ebs << std::endl;
  std::cout << "Erep = " << erep << std::endl;
  std::cout << "Etot = " << etot << std::endl;

  return 0;
}

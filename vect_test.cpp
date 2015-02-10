#include <iostream>
#include <vector>
#include <cmath>
#include "hamiltonian.h"
#include "functions.h"

int main() {
  int n=3;
  std::vector<int> type(n);
  std::vector<double> posx(n);
  std::vector<double> posy(n);
  std::vector<double> posz(n);
  int i,j;
  std::vector<double> eigvects(16*n*n);

  type.at(0) = 6;
  type.at(1) = 6;
  type.at(2) = 6;
  posx.at(0) = -1.4000;
  posx.at(1) = 0.0000;
  posx.at(2) = 1.6000;
  posy.at(0) = 0.0000;
  posy.at(1) = 0.3000;
  posy.at(2) = 0.0000;
  posz.at(0) = 0.0000;
  posz.at(1) = 0.0000;
  posz.at(2) = 0.0000;

  double ebs = Hamiltonian(n, &type, &posx, &posy, &posz, &eigvects);
  double erep = Erep(type, posx, posy, posz);  
  double etot = ebs + erep;

  std::cout << "Ebs = " << ebs << std::endl;
  std::cout << "Erep = " << erep << std::endl;
  std::cout << "Etot = " << etot << std::endl;

  return 0;
}

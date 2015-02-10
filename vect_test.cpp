#include <iostream>
#include <vector>
#include <cmath>
#include "hamiltonian.h"
#include "functions.h"

int main() {
  int n=2;
//  std::vector<int> type(n);
  std::vector<double> posx(n);
  std::vector<double> posy(n);
  std::vector<double> posz(n);
  int i,j;
  std::vector<double> eigvects(16*n*n);

//  type.at(0) = 6;
//  type.at(1) = 6;
//  type.at(2) = 6;
  posx.at(0) = 0.0000;
  posx.at(1) = 1.3000;
//  posx.at(2) = 2.4800;
  posy.at(0) = 0.0000;
  posy.at(1) = 0.0000;
//  posy.at(2) = 0.0000;
  posz.at(0) = 0.0000;
  posz.at(1) = 0.0000;
//  posz.at(2) = 0.0000;

  double ebs = Hamiltonian(n,&posx, &posy, &posz, &eigvects);
  double erep = Erep(&posx, &posy, &posz);  
  double etot = ebs + erep;

  std::cout << "Ebs = " << ebs << std::endl;
  std::cout << "Erep = " << erep << std::endl;
  std::cout << "Etot = " << etot << std::endl;

  return 0;
}

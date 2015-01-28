#include <iostream>
#include <vector>
#include <cmath>
#include "hamiltonian.h"

int main() {
  std::vector<int> type(2);
  std::vector<double> posx(2);
  std::vector<double> posy(2);
  std::vector<double> posz(2);
  int n = 2;
  type.at(0) = 6;
  type.at(1) = 6;
  posx.at(0) = 0.0000;
  posx.at(1) = 1.400;
  posy.at(0) = 0.0000;
  posy.at(1) = 0.0000;
  posz.at(0) = 0.0000;
  posz.at(1) = 0.0000;

  Hamiltonian(n, &type, &posx, &posy, &posz);

  return 0;
}

#include <iostream>
#include <vector>
#include <cmath>
#include "hamiltonian.h"

/*void Hamiltonian(int n, std::vector<int>* type, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz);*/

int main() {
  std::vector<int> type(2);
  std::vector<double> posx(2);
  std::vector<double> posy(2);
  std::vector<double> posz(2);
  int n = 2;
  type.at(0) = 6;
  type.at(1) = 6;
  posx.at(0) = 0.0000;
  posx.at(1) = 1.2400;
  posy.at(0) = 0.0000;
  posy.at(1) = 0.0000;
  posz.at(0) = 0.0000;
  posz.at(1) = 0.0000;

  std::cout << "Checkpoint 0" << std::endl;

  Hamiltonian(n, &type, &posx, &posy, &posz);
  //Test(n, &posx);

  /*H_MD.at(0) = 0.00382;
  H_MD.at(1) = 2.2398;
  H_MD.at(2) = 4.00382;
  double var = H_MD[0];
  std::cout << "var = " << var << std::endl;
  double var2 = H_MD.at(2);
  std::cout << "var2 = " << var2 << std::endl;*/
  return 0;
}

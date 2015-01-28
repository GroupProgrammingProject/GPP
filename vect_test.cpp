#include <iostream>
#include <vector>
#include <cmath>

int main() {
  std::vector<double> H_MD(3);
 H_MD.at(0) = 0.00382;
 H_MD.at(1) = 2.2398;
 H_MD.at(2) = 4.00382;
 double var = H_MD[0];
 std::cout << "var = " << var << std::endl;
 double var2 = H_MD.at(2);
 std::cout << "var2 = " << var2 << std::endl;
}

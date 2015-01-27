#include <cmath>
#include <iostream>
#include <iomanip>
#include "functions.h"

using namespace std;

int main() {
 
  hello();
  double sos=s(2.45);
  cout <<"value cutoff s "<<sos<<endl;
  sos=s(2.449999999);
  cout <<"value cutoff s "<<sos<<endl;

  double ooo=o(2.57);
  cout <<"value cutoff o "<<ooo<<endl;
  ooo=o(2.569999999999);
  cout <<"value cutoff o "<<ooo<<endl;


};

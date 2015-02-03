#include <cmath>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include "functions.h"

using namespace std;

int main() {
 
  hello();
  double sos;

  /*ofstream tm;
  tm.open ("scaling_s(r).txt");

  for(double i=2.4;i<2.6;i=i+0.001){
  sos=s(i,6,6);
  tm<<i<<"\t"<<sos<<endl;
  };
  tm.close();
  */
  sos=s(2.449999999,6,6);
  cout <<"value cutoff s "<<sos<<endl;

  double ooo=o(2.57);
  cout <<"value cutoff o "<<ooo<<endl;
  ooo=o(2.569999999999);
  cout <<"value cutoff o "<<ooo<<endl;


 std:vector<double> wekt;
  for(int i=0;i<10;i++){
    wekt.push_back(3.32+i*0.04);
  }
  test(wekt);

};

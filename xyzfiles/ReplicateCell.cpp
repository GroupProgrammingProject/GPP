#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "../include/readinxyz.h"
#include "../include/ScaleGeom.h"

int main(int argc, char* argv[]){
  if (argc<2){std::cout<<"You should append two files to the main object!"<<std::endl;}
  if (argc!=3){std::cout<<"You should append one input and one output xyz files to the main!!"<<std::endl;}
  
  std::vector<double> posx, posy, posz;
  ReadInXYZ (argv[1], &posx, &posy, &posz);
  // Number of atoms
  int n=posx.size();

  // Replicate cell in n dimensions
  double ascale = 4;
  double bscale = 1;
  double cscale = 1;
  int n_new = ascale*bscale*cscale*n;
  double a = 20;
  double b = 20;               
  double c = 7.37854;
  std::vector<double> posxnew(n_new), posynew(n_new), posznew(n_new); 
  ScaleCell(n, ascale, bscale, cscale, a, b, c, &posx, &posy, &posz, &posxnew, &posynew, &posznew);
  
  // Print output .xyz cell
  std::string outputfile = argv[2];
  std::ofstream cell(outputfile.c_str());
  cell << n_new << "\n";
  cell << outputfile << "\t" << ascale*a << "\t" << bscale*b << "\t" << cscale*c << "\n";
  for (int j=0;j<n_new;j++){
	 cell << "C\t" << std::setprecision(6) << posxnew.at(j) << "\t" << std::setprecision(6) << posynew.at(j) <<  "\t" << std::setprecision(6) << posznew.at(j) << "\n";
  }

  return 0;
}

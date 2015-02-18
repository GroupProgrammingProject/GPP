#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "readinxyz.h"
#include "ScaleGeom.h"

int main(int argc, char* argv[]){
  if (argc<2){std::cout<<"You should append two files to the main object!"<<std::endl;}
  if (argc!=3){std::cout<<"You should append one input and one output xyz files to the main!!"<<std::endl;}
  
  std::vector<double> posx, posy, posz;
  ReadInXYZ (argv[1], &posx, &posy, &posz);
  // Number of atoms
  int n=posx.size();

  // Replicate cell in n dimensions
  double ascale = 0.9483940421888254;
  double bscale = 0.9483940421888254;
  double cscale = 0.9483940421888254;
  int n_new = n;
  std::vector<double> posxnew(n_new), posynew(n_new), posznew(n_new); 
  ScaleGeom(n, ascale, bscale, cscale, &posx, &posy, &posz, &posxnew, &posynew, &posznew);
  
  // Print output .xyz cell
  std::string outputfile = argv[2];
  std::ofstream cell(outputfile.c_str());
  cell << n_new << "\n";
  cell << outputfile << "\n";
  for (int j=0;j<n_new;j++){
	 cell << "C\t" << std::setprecision(6) << posxnew.at(j) << "\t" << std::setprecision(6) << posynew.at(j) <<  "\t" << std::setprecision(6) << posznew.at(j) << "\n";
  }

  return 0;
}

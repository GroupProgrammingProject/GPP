#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include "../include/readinxyz.h"
#include "../include/ScaleGeom.h"

// Takes "input.xyz output.xyz A B C" as input where A, B and C are the integers by which you wish to replicate the input.xyz structure

int main(int argc, char* argv[]){
  if (argc<5){std::cout<<"You should append two files and the numbers by which you wish to replicate the cell along each axis to the main object!"<<std::endl;}
  if (argc!=6){std::cout<<"You should append one input and one output xyz files and the numbers by which you wish to replicate the cell along each axis to the main!!"<<std::endl;}
  
  bool pbc = 1;
  std::vector<bool> velspec;
  std::vector<double> lats(3);
  std::vector<double> posx, posy, posz, vxin, vyin, vzin;
  ReadInXYZ (argv[1], &posx, &posy, &posz, &vxin, &vyin, &vzin, &lats, pbc, &velspec);
  // Number of atoms
  int n=posx.size();

  // Replicate cell in n dimensions
  int scale[3] = {atoi(argv[3]),atoi(argv[4]),atoi(argv[5])};
  std::cout << "scale: " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
  int n_new = scale[0]*scale[1]*scale[2]*n;
  std::cout << "n_new = " << n_new << std::endl;
  std::vector<double> posxnew(n_new), posynew(n_new), posznew(n_new), latsnew(3); 
  ScaleCell(n, scale, &lats, &posx, &posy, &posz, &posxnew, &posynew, &posznew, &latsnew, pbc);

  // Print output .xyz cell
  std::string outputfile = argv[2];
  std::ofstream cell(outputfile.c_str());
  cell << n_new << "\n";
  cell << outputfile << "\t" << latsnew[0] << "\t" << latsnew[1] << "\t" << latsnew[2] << "\n";
  for (int j=0;j<n_new;j++){
	 cell << "C\t" << std::setprecision(6) << posxnew.at(j) << "\t" << std::setprecision(6) << posynew.at(j) <<  "\t" << std::setprecision(6) << posznew.at(j) << "\n";
  }

  return 0;
}

#include <iostream>
#include <fstream>
#include <vector>
#include "readinxyz.h"
#include "ScaleGeom.h"

int main(int argc, char* argv[]){
  if (argc<1){std::cout<<"You should append a file to the main object!"<<std::endl;}
  if (argc!=2){std::cout<<"You should append one and only one xyz file to the main!!"<<std::endl;}
  
  std::vector<int> type;
  std::vector<double> posx, posy, posz;
  ReadInXYZ (argv[1], &type, &posx, &posy, &posz);
  // Number of atoms
  int n=posx.size();

  /*// Scale cell size
  int n_new = n;
  double ascale = 1.2;                                                  // By what float value to scale the x axis
  double bscale = 1.2;                                                  // By what float value to scale the y axis
  double cscale = 1.2;                                                  // By what float value to scale the z axis
  std::vector<double> posxnew(n), posynew(n), posznew(n);               // Define new vectors to hold new atom positions
  ScaleGeom(n, ascale, bscale, cscale, &posx, &posy, &posz, &posxnew, &posynew, &posznew);*/

  // Replicate cell in n dimensions
  double ascale = 4;                                                    // By how many times to replicate cell in x dimension
  double bscale = 4;                                                    // By how many times to replicate cell in y dimension
  double cscale = 4;                                                    // By how many times to replicate cell in z dimension
  int n_new = ascale*bscale*cscale*n;                                   // Number of atoms in new unit cell
  double a = 3.57*0.5;                                                  // a lattice parameter (in angstroms)
  double b = 3.57*0.5;                                                  // b lattice parameter (in angstroms)
  double c = 3.57*0.5;                                                  // c lattice parameter (in angstroms)
  std::vector<double> posxnew(n_new), posynew(n_new), posznew(n_new);   // Define new vectors to hold new atom positions
  ScaleCell(n, ascale, bscale, cscale, a, b, c, &posx, &posy, &posz, &posxnew, &posynew, &posznew);
  
  // Print output .xyz cell
  std::ofstream cell("diamondcell_out.xyz");
  cell << n_new << "\n";
  cell << "Ouput name" << "\n";
  for (int j=0;j<n_new;j++){
	 cell << "C\t" << posxnew.at(j) << "\t" << posynew.at(j) <<  "\t" << posznew.at(j) << "\n";
  }

  return 0;
}

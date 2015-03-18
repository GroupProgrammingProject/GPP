#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include "../include/readinxyz.h"
#include "../include/kpointsfunctions.h"

// Main for generating an orthorhombic grid of kpoints
// Takes "input.xyz output.kpts A B C" where input.xyz is a structure file, output.kpts is a file to hold a set of points in kspace to be read by the calculation main and A, B and C are the number of kpoints to be generated along the reciprocal a, b and c directions

int main(int argc, char* argv[]) {
  if (argc<5){std::cout<<"You should append two files and 3 kpoint grids to the main object!"<<std::endl;}
  if (argc!=6){std::cout<<"You should append one xyz, one .kpts file and a kpoint grid to the main!!"<<std::endl;}
  std::vector<double> posx, posy, posz;                    // Required for ReadInXYZ
  std::vector<double> lats(3);                             // Vector to take lattice parameters
  bool pbc = 1;                                            // If generating a kpoint grid then using PBCs
  ReadInXYZ (argv[1], &posx, &posy, &posz, &lats, pbc);    // Read in lattice parameters from .xyz file
  int kn[3] = {atoi(argv[3]),atoi(argv[4]),atoi(argv[5])}; // kpoint grid read in from command line
  genkgrid(argv[2],&lats,kn,0);                            // Generate grid in specified kpoint file
}

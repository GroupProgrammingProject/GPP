#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include "../include/readinxyz.h"
#include "../include/kpointfunctions.h"

// Main for generating a path between two points in kspace
// Takes "input.xyz output.kpts P0 Q0 R0 P1 Q1 R1 npts" where input.xyz is a structure file, output.kpts is a file to hold a set of points in kspace to be read by the calculation main and P, Q and R are doubles representing coordinates in kspace (along reciprocal a, b and c directions) for the inital and final points on the staight path. npts is then the number of points to be generated along this path.

int main(int argc, char* argv[]) {
  if (argc<9){std::cout<<"You should append two files an initial and final kpoint and the number of points to be generated to the main object!"<<std::endl;}
  if (argc!=10){std::cout<<"You should append one xyz, one .kpts file, the coordinates of initial and final kpoints and the number of gridpoints to the main!!"<<std::endl;}
  std::vector<double> posx, posy, posz;                         // Required for ReadInXYZ
  std::vector<double> lats(3);                                  // Vector to take lattice parameters
  bool pbc = 1;                                                 // If generating a kpoint grid then using PBCs
  ReadInXYZ (argv[1], &posx, &posy, &posz, &lats, pbc);         // Read in lattice parameters from .xyz file
  double kpt0[3] = {atof(argv[3]),atof(argv[4]),atof(argv[5])}; // Inital kpoint read in from command line 
  double kpt1[3] = {atof(argv[6]),atof(argv[7]),atof(argv[8])}; // Final kpoint read in from command line 
  int npts = atoi(argv[9]);                                     // Number of points along path
  genkpath(argv[2], &lats, kpt0, kpt1, npts);                   // Generate path in specified kpoint file
}

#ifndef SCALEGEOM_H
#define SCALEGEOM_H

#include <cmath>
#include <iostream>
#include <vector>

void ScaleGeom(int n, double scale[3], std::vector<double>* lats, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew, std::vector<double>* latsnew, bool pbc);
void ScaleCell(int n, int scale[3], std::vector<double>* lats, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew, std::vector<double>* latsnew, bool pbc);

// Scale the positions of atoms along each of three axes
void ScaleGeom(int n, double scale[3], std::vector<double>* lats, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew, std::vector<double>* latsnew, bool pbc) {
  for (int i=0;i<n;i++) {                                          // For each atom in original cell
	 (*posxnew).at(i) = scale[0]*(*posx).at(i);                     // Scale in each dimension
	 (*posynew).at(i) = scale[1]*(*posy).at(i);
	 (*posznew).at(i) = scale[2]*(*posz).at(i);
  }
  // Scale lattice parameters if using PBCs
  if (pbc == 1) {
	 for (int i=0;i<3;i++) {
		(*latsnew).at(i) = scale[i]*(*lats).at(i);
	 }
  }
}

// Replicate the existing cell an integer number of times along three axes
void ScaleCell(int n, int scale[3], std::vector<double>* lats, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew, std::vector<double>* latsnew, bool pbc) {
  int j = 0;
  for (int i=0;i<n;i++) {                                          // For each atom in original cell
	 for (int jx=0;jx<scale[0];jx++) {                              // Replicate atom in each dimension
		for (int jy=0;jy<scale[1];jy++) {
		  for (int jz=0;jz<scale[2];jz++) {
			 (*posxnew).at(j) = jx*(*lats).at(0) + (*posx).at(i);
			 (*posynew).at(j) = jy*(*lats).at(1) + (*posy).at(i);
			 (*posznew).at(j) = jz*(*lats).at(2) + (*posz).at(i);
			 std::cout << "New atom " << j << " at " << (*posxnew).at(j) << " " << (*posynew).at(j) << " " << (*posznew).at(j) << std::endl;
			 j++;
		  }
		}
	 }
  }
  // Scale lattice parameters if using PBCs
  if (pbc == 1) {
	 for (int i=0;i<3;i++) {
		(*latsnew).at(i) = scale[i]*(*lats).at(i);
	 }
  }
}

#endif

#ifndef SCALEGEOM_H
#define SCALEGEOM_H

#include <cmath>
#include <iostream>
#include <vector>

void ScaleGeom(int n, double ascale, double bscale, double cscale, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew);
void ScaleCell(int n, int ascale, int bscale, int cscale, double a, double b, double c, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew);

// Scale the positions of atoms along each of three axes
void ScaleGeom(int n, double ascale, double bscale, double cscale, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew) {
  for (int i=0;i<n;i++) {
	 (*posxnew).at(i) = ascale*(*posx).at(i);
	 (*posynew).at(i) = bscale*(*posy).at(i);
	 (*posznew).at(i) = cscale*(*posz).at(i);
  }
}

// Replicate the existing cell an integer number of times along three axes
void ScaleCell(int n, int ascale, int bscale, int cscale, double a, double b, double c, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew) {
  int j;
  for (int i=0;i<n;i++) {
	 for (int jx=0;jx<ascale;jx++) {
		for (int jy=0;jy<bscale;jy++) {
		  for (int jz=0;jz<cscale;jz++) {
			 (*posxnew).at(j) = jx*a + (*posx).at(i);
			 (*posynew).at(j) = jy*b + (*posy).at(i);
			 (*posznew).at(j) = jz*c + (*posz).at(i);
			 std::cout << "New atom " << j << " at " << (*posxnew).at(j) << " " << (*posynew).at(j) << " " << (*posznew).at(j) << std::endl;
			 j++;
		  }
		}
	 }
  }
}

#endif

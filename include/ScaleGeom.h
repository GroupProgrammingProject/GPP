#ifndef SCALEGEOM_H
#define SCALEGEOM_H

#include <cmath>
#include <iostream>
#include <vector>

void ScaleGeom(int n, double scale[3], std::vector<double>* lats, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew, std::vector<double>* latsnew, bool pbc);
void ScaleCell(int n, int scale[3], std::vector<double>* lats, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew, std::vector<double>* latsnew, bool pbc);

#endif

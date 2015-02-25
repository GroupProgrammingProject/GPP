#ifndef SCALEGEOM_H
#define SCALEGEOM_H

#include <cmath>
#include <iostream>
#include <vector>

void ScaleGeom(int n, double ascale, double bscale, double cscale, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew);
void ScaleCell(int n, int ascale, int bscale, int cscale, double a, double b, double c, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* posxnew, std::vector<double>* posynew, std::vector<double>* posznew);

#endif

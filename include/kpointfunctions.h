#ifndef KPOINTFUNCTIONS_H
#define KPOINTFUNCTIONS_H

#include <vector>
#include <cmath>
#include <fstream>

void readinkpoints(char* filename, std::vector<std::vector<double> >* kpoints);
void genkgrid(char* filename, std::vector<double>* lats, int kgrid[3], bool gamma);
void genkpath(char* filename, std::vector<double>* lats, double kpt0[3], double kpt1[3], int npts);

#endif

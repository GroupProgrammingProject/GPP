#ifndef KPOINTFUNCTIONS_H
#define KPOINTFUNCTIONS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include "functions.h"
#include "MolDyn.h"

void readinkpoints(char* filename, std::vector<std::vector<double> >* kpoints);
void genkgrid(char* filename, std::vector<double>* lats, int kgrid[3], bool gamma);
void genkpath(char* filename, std::vector<double>* lats, double kpt0[3], double kpt1[3], int npts);
void kforces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXcd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<std::complex<double> >* fx, std::vector<std::complex<double> >* fy, std::vector<std::complex<double> >* fz, std::vector<double>* kvec);
double kHamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum);

#endif

#ifndef KPOINTFUNCTIONS_H
#define KPOINTFUNCTIONS_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <time.h>
#include <stdlib.h>
#include <cstdlib> 
#include <stdio.h>
#include <vector>
#include <ctime>
#include "band_hamiltonian.h"
#include "geometryinfo.h"
#include "functions.h"
//#include "MolDyn.h"
#include "nr3.h"
#include "ran.h"

void readinkpoints(char* filename, std::vector<std::pair<std::vector<double>,double> >* kpoints, bool ksymm);
void genkgrid(char* filename, std::vector<double>* lats, int kgrid[3], bool gamma);
void genkpath(char* filename, std::vector<double>* lats, double kpt0[3], double kpt1[3], int npts);
void kforces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXcd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz, std::vector<double>* kvec, std::vector<double>* TBparam);
double kHamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum, std::vector<double>* TBparam);
double avekenergy(int N,int norbs,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr,std::vector<std::pair<std::vector<double>,double> >* kpoints,std::vector<double>* TBparam);
double avekforces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz, std::vector<std::pair<std::vector<double>,double> >* kpoints,std::vector<double>* TBparam);
double kverlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, double &ebs, std::vector<double>* lats, bool pbc,double T, double nu,bool ander, std::vector<std::pair<std::vector<double>,double> >* kpoints,std::vector<double>* TBparam);

#endif

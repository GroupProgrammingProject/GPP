#ifndef PHONONS_H
#define PHONONS_H

#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "MolDyn.h"
#include "hamiltonian.h"

void normalmodes(int n,int norbs,double rc,double m,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* eigvects, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz, std::vector<double>* eigfreq,std::vector<double>* TBparam);
/*Inputs, in order: #atoms; #orbitals; cut-off radius; mass of atoms;atom vector distances; modulus of vector distances; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; forces vectors,eigenmodes*/

#endif

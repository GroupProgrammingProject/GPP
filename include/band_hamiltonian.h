#ifndef BAND_HAMILTONIAN_H
#define BAND_HAMILTONIAN_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <cmath>
#include <complex>
#include "Gethijab.h"
#include "functions.h"

double band_Hamiltonian(int n,Eigen::MatrixXd* modr, Eigen::MatrixXd* rx,Eigen::MatrixXd* ry,Eigen::MatrixXd* rz,Eigen::MatrixXcd* eigvects,std::vector<double>* eigenvalaar, std::vector<double>* k, bool verbose);
double dot(std::vector<double>* a, std::vector<double>* b);

#endif

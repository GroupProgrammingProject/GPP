#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "Gethijab.h"
#include "functions.h"

double Hamiltonian(int n,Eigen::MatrixXd* modr, Eigen::MatrixXd* rx,Eigen::MatrixXd* ry,Eigen::MatrixXd* rz,Eigen::MatrixXd* eigvects, bool verbose);
#endif

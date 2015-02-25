#ifndef COMPLEX_SOLVER_HAMN_H
#define COMPLEX_SOLVER_HAMN_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "Gethijab.h"
#include "functions.h"

double Hamiltonian(int n, std::vector<int>* type, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* eigvects);

#endif

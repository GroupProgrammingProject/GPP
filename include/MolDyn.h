#ifndef MOLDYN_H
#define MOLDYN_H

// Functionality for molecular dynamics
#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <cstdlib> 
#include <stdio.h>
#include <vector>
#include <ctime>
#include "functions.h"
#include "Gethijab.h"
#include "geometryinfo.h"
#include <vector>

//Vector syntax: std::vector<double>. Assignments: xi=(*v).at(i). Inputs as pointer: std::vector<double>*

/*Change vector syntax as: std::vector<double>. Assignments xi=(*v).at(i). Inputs as pointer: std::vector<double>*, to call in functs ad &*/

double verlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>*vy, std::vector<double>* vz, Eigen::MatrixXd* c,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr);/*Inputs, in order: #orbitals; timestep; xyz arrays; velocity arrays; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; atom vector distances; modulus of vector distances.*/

void forces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz);/*Inputs, in order: #atoms; #orbitals; atom vector distances; modulus of vector distances; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; forces vectors.*/

void velocity(double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T);/*Inputs, in order: #atoms; mass; vectors of velocities; temperature.*/

double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum);/*Inputs, in order: first atom index; second atom index; first orbital index; second orbital index; vector of direction cosines; modulus of distance between i and j; switch for component.*/ 
#endif

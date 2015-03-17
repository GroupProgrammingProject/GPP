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
#include <Eigen/Dense>
#include "functions.h"
#include "hamiltonian.h"
#include "Gethijab.h"
#include "geometryinfo.h"
#include "nr3.h"
#include "ran.h"

//Vector syntax: std::vector<double>. Assignments: xi=(*v).at(i). Inputs as pointer: std::vector<double>*

/*Change vector syntax as: std::vector<double>. Assignments xi=(*v).at(i). Inputs as pointer: std::vector<double>*, to call in functs ad &*/

double verlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>*vy, std::vector<double>* vz, Eigen::MatrixXd* c,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr,double &ebs, std::vector<double>* lats, bool pbc,double T,double nu,bool ander);
/*Inputs, in order: #orbitals; timestep; xyz arrays; velocity arrays; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; atom vector distances; modulus of vector distances.*/

void forces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz);/*Inputs, in order: #atoms; #orbitals; atom vector distances; modulus of vector distances; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; forces vectors.*/

void velocity(double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T);/*Inputs, in order: #atoms; mass; vectors of velocities; temperature.*/

double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum);/*Inputs, in order: first atom index; second atom index; first orbital index; second orbital index; vector of direction cosines; modulus of distance between i and j; switch for component.*/ 

int GeomOpt(int norbs,double rc,double rv,double m,double dt,int nmd, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz,std::vector<double>* refposx, std::vector<double>* refposy, std::vector<double>* refposz, Eigen::MatrixXd* eigvects,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, std::vector<double>* lats, bool pbc,double T,double nu,double h,bool verb,int nprint);
/*Inputs, in order: #orbitals; cut-off radius; Verlet radius; atomic mass; timestep; number of simulated annealing steps; xyz arrays; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; atom vector distances; modulus of vector distances; boundary conditions vector; pbc switch; initial temperature; Andersen frequency; steepest descent step; output switch.*/

//void chain(int N,double m,double T,double dt,double q1,double q2,double &xi1,double &xi2,double &vxi1,double &vxi2,std::vector<double>* vx,std::vector<double>* vy,std::vector<double>* vz,double &kin);
#endif

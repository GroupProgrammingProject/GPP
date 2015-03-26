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
#include "hamiltonian.h"
#include "functions.h"
#include "Gethijab.h"
#include "geometryinfo.h"
#include <vector>
#include "nr3.h"
#include "ran.h"
#include "kpointsfunctions.h"

/*verlet performs a single MD step using Velocity Verlet integration. It's possible to decide whether to use PBCs and/or Andersen's thermostat.*/ 
double verlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>*vy, std::vector<double>* vz, Eigen::MatrixXd* c,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr,double &ebs, std::vector<double>* lats, bool pbc,double T,double nu,bool ander,std::vector<double>* TBparam);
/*Inputs, in order: #orbitals (norbs); cut-off radius (rc); Verlet radius (rv); atomic mass (m); timestep [ps] (dt); position arrays (x,y,z); reference position arrays (refx,refy,refz); velocity arrays (vx,vy,vz); matrix of eigenvectors, with vectors as rows (c); nearest neighbour lists (vector nnear and matrix inear); matrices of vector distances (rx,ry,rz); matrix of moduli of vector distances (modr); band structure energy (ebs); supercell dimensions (lats); periodic boundary conditions switch (pbc); temperature (T); frequency of velocity rescalings for thermostat (nu); thermostat switch (ander); tight binding parameters (TBparam).*/

/*forces computes both band structure and repulsive forces analytically*/
void forces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz,std::vector<double>* TBparam);
/*Inputs, in order: #atoms; #orbitals; atom vector distances; modulus of vector distances; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; forces vectors.*/

/*velocity assigns initial velocities to all the atoms according to the initial temperature.*/
void velocity(double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T);/*Inputs, in order: mass; vectors of velocities; temperature.*/

/*Hamder returns analytically detivatives of the (unscaled) matrix elements of the TB matrix.*/
double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum, std::vector<double>* TBparam);
/*Inputs, in order: first atom index; second atom index; first orbital index; second orbital index; vector of direction cosines; modulus of distance between i and j; switch for component.*/ 

/*GeomOpt performs gemoetrical optimisation of a given structure using simulated annealing followed by steepest descent.*/
int GeomOpt(int norbs,double rc,double rv,double m,double dt,int nmd, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz,std::vector<double>* refposx, std::vector<double>* refposy, std::vector<double>* refposz, Eigen::MatrixXd* eigvects,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, std::vector<double>* lats, bool pbc,double T,double nu,double h,bool verb,int nprint,std::vector<double>* TBparam, double tol,int maxsteep,bool kpts,std::vector<std::pair<std::vector<double>,double> >* kpoints);
/*Inputs, in order: #orbitals; cut-off radius; Verlet radius; atomic mass; timestep; number of simulated annealing steps; xyz arrays; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; atom vector distances; modulus of vector distances; boundary conditions vector; pbc switch; initial temperature; Andersen frequency; steepest descent step; output switch.*/

#endif

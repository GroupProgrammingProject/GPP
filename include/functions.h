#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

//scaling functions, notations as in the paper 
double ts (double r);
double dts (double r,double dx);
double to (double r);
double dto (double r, double dx);
double s (double r);
double ds (double r,double dx);
double o (double r);
double f0(double x);
double d_f0(double x);
double d_o(double r,double dx);
double X (Eigen::MatrixXd* modr,int n, int i );
double Erep (Eigen::MatrixXd* modr);

#endif

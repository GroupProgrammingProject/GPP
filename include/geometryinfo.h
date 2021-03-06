#ifndef GEOMETRYINFO_H
#define GEOMETRYINFO_H

#include <vector>
#include <cmath>
#include "vectorfunctions.h"
#include <Eigen/Dense>

// Function declarations
double Abs(double c);
void DistanceComp(Eigen::MatrixXd* r, std::vector<double>* pos);
void GetDistanceComponents(Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz);
void GetDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* lats, double rv, bool pbc);
void GetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz);
void GetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc);
void PbcDistanceComp(Eigen::MatrixXd* r, std::vector<double>* pos, double latticeconst, double rv);
void PbcGetDistanceComponents(Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double a, double b, double c, double rv);
void PbcGetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* lats, double rv);
void NearestNeighbours(Eigen::MatrixXi* inear, std::vector<int>* nnear, Eigen::MatrixXd* modr, double rv);
bool RecalculateNearestNeighbours(std::vector<double>* refposx, std::vector<double>* refposy, std::vector<double>* refposz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc, double rv);
void SetEqual(std::vector<double>* a, std::vector<double>* b);
void pbcshift(std::vector<double>* xvect,std::vector<double>* yvect,std::vector<double>* zvect,std::vector<double>* lats);

#endif

#include <cmath>

void ScaleGeom(int n, double ascale, double bscale, double cscale, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz) {
  for (int i=0;i<n;i++) {
	 (*posx).at(i) = ascale*(*posx).at(i);
	 (*posy).at(i) = bscale*(*posy).at(i);
	 (*posz).at(i) = cscale*(*posz).at(i);
  }
}

void ScaleCell(int n, int ascale, int bscale, int cscale, double a, double b, double c, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz) {
  for (int i=0;i<n;i++) {
	 for (int j=0;j<ascale;j++) {
		(*posx).at(i+n*j) = j*a + (*posx).at(i);}
	 for (int j=0;j<bscale;j++) {
		(*posy).at(i+n*j) = j*b + (*posy).at(i);}
	 for (int j=0;j<cscale;j++) {
		(*posz).at(i+n*j) = j*c + (*posz).at(i);}
  }
}

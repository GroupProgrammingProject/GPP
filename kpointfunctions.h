#ifndef KPOINTFUNCTIONS_H
#define KPOINTFUNCTIONS_H

#include <vector>
#include <cmath>
#include <fstream>

void readinkpoints(char* filename, std::vector<std::vector<double> >* kpoints) {
  std::ifstream infile(filename);
  double k0, k1, k2;
  std::vector<double> kvec(3);
  while (infile>>k0>>k1>>k2){
	 kvec.at(0) = k0;
	 kvec.at(1) = k1;
	 kvec.at(2) = k2;
	 (*kpoints).push_back(kvec);
  }
}

void generatekpoints(char* filename, std::vector<double>* lats, int kgrid[3], bool gamma) {
  std::ofstream outfile(filename);
  std::cout << "kpoint grid " << kgrid[0] << " " << kgrid[1] << " " << kgrid[2] << std::endl;
  double kstep[3] = {2*M_PI/((*lats).at(0)*kgrid[0]), 2*M_PI/((*lats).at(1)*kgrid[1]), 2*M_PI/((*lats).at(2)*kgrid[2])}; 
  double kvec[3];
  for(double k0=-M_PI/(*lats).at(0);k0<M_PI/(*lats).at(0)-0.5*kstep[0];k0=k0+kstep[0]){
	 for(double k1=-M_PI/(*lats).at(1);k1<M_PI/(*lats).at(1)-0.5*kstep[1];k1=k1+kstep[1]){
		for(double k2=-M_PI/(*lats).at(2);k2<M_PI/(*lats).at(2)-0.5*kstep[2];k2=k2+kstep[2]){
		  kvec[0] = k0 + gamma*0.5*kstep[0];
		  kvec[1] = k1 + gamma*0.5*kstep[1];
		  kvec[2] = k2 + gamma*0.5*kstep[2];
		  outfile << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << "\n";
		}
	 }
  }
}

#endif

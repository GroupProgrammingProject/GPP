#include "../include/kpointsfunctions.h"
// Read in kpoints from a file (e.g. for a calculation)
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

// Generate a user specified grid of kpoints (only for orthorhombic symmetries)
void genkgrid(char* filename, std::vector<double>* lats, int kgrid[3], bool gamma) {
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

// Generate a path of npts kpoints between two user specified points in kspace
void genkpath(char* filename, std::vector<double>* lats, double kpt0[3], double kpt1[3], int npts) {
  std::ofstream outfile(filename);
  double kstep[3];
  double kvec[3];
  for (int i=0;i<3;i++) {                                        // For each direction
	 kstep[i] = (kpt1[i] - kpt0[i])/(double)npts;                 // kstep is difference divided by npts
  }
  for (jnt j=0;j<npts;j++) {
	 for (int i=0;i<3;i++) {
		kvec[i] = kpt0[i] + j*kstep[i];                            // Iterate along steps for desired number of kpoints
	 }
	 outfile << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << "\n";
  }
}



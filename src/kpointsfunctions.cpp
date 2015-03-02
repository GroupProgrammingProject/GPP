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
  for (int j=0;j<npts;j++) {
	 for (int i=0;i<3;i++) {
		kvec[i] = kpt0[i] + j*kstep[i];                            // Iterate along steps for desired number of kpoints
	 }
	 outfile << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << "\n";
  }
}

void kforces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz)
{ int k,i,j,l,lp,n,m,nearlabel; /* dummy indeces for cycles*/
  std::vector<double> ddnorm(3);
  double sumphinn,sumphi,derivx,derivy,derivz,dx,dy,dz,r;
  double dualeigen;

  for(i=0;i<N;i++){ /*initialisation of forces*/
    (*fx).at(i)=0;
    (*fy).at(i)=0;
    (*fz).at(i)=0;
  }
  
  for(i=0;i<N;i++){ /*Cycle to compute band structure forces on atom i*/
    sumphi=0;
    for(k=0;k<(*nnear).at(i);k++){
      nearlabel=(*inear)(i,k);
      r=(*modr)(i,nearlabel);
      if(r<rc){
	sumphi=sumphi+o(r);
      }
    }
    for(j=0;j<(*nnear).at(i);j++){ /*Cycle spanning the nearest neighbours of i*/
      nearlabel=(*inear)(i,j);
      sumphinn=0;
      for(m=0;m<(*nnear).at(nearlabel);m++){     
	dx=(*rx)(nearlabel,(*inear)(nearlabel,m)); /*Definition of vector distances*/
	dy=(*ry)(nearlabel,(*inear)(nearlabel,m));
	dz=(*rz)(nearlabel,(*inear)(nearlabel,m));
	r=(*modr)(nearlabel,(*inear)(nearlabel,m)); /*Modulus of distance*/
	if(r<rc){
	  sumphinn=sumphinn+o(r);
	}
      }
      dx=(*rx)(i,nearlabel); /*Definition of vector distances*/
      dy=(*ry)(i,nearlabel);
      dz=(*rz)(i,nearlabel);
      
      r=(*modr)(i,nearlabel); /*Modulus of distance*/
      if(r<rc){
	ddnorm.at(0)=dx/r;
	ddnorm.at(1)=dy/r;
	ddnorm.at(2)=dz/r;
      
	for(l=0;l<norbs;l++){ /*Cycle spanning the first orbital type*/
	  for(lp=0;lp<norbs;lp++){ /*Cycle spanning the second orbital type*/
	    derivx=ds(r,dx)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*kHamder(i,nearlabel,l,lp,&ddnorm,r,0);
	    derivy=ds(r,dy)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*kHamder(i,nearlabel,l,lp,&ddnorm,r,1);
	    derivz=ds(r,dz)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*kHamder(i,nearlabel,l,lp,&ddnorm,r,2);
	    for(n=0;n<norbs*N;n++){ /*Cycle spanning the level of the eigenvector*/
	      dualeigen=(*c)(n,l+i*norbs)*(*c)(n,lp+nearlabel*norbs);  
	      (*fx).at(i)=(*fx).at(i)-2*derivx*dualeigen;
	      (*fy).at(i)=(*fy).at(i)-2*derivy*dualeigen;
	      (*fz).at(i)=(*fz).at(i)-2*derivz*dualeigen;
	    }
	  }
	}
	//calculation of repulsive forces
	(*fx).at(i)=(*fx).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dx);
	(*fy).at(i)=(*fy).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dy);
	(*fz).at(i)=(*fz).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dz);
      }
    }
  }  
}

//kHamder() returns value of Hamilatonian matrix element differentiated wrt x,y or z.
double kHamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum){
  int k; //for looping
  double h,V[4];//h,Es,Ep and V[4] is only used locally in Gethijab_der()
  V[0]=-5;V[1]=4.7;V[2]=5.5;V[3]=-1.55;//CC interaction 0=ss_sigma, 1=sp_sigma, 2=pp_sigma, 3=pp_pi
  
  int vconum[3]; //array to hold label of x,y,z. Only the coord we are differentiating wrt is non-zero.
  vconum[0]=0;
  vconum[1]=0;
  vconum[2]=0;
  vconum[conum]=1;

  //start V&G routine
  if(i==j){h=0;}
  else if(a*b==0){
    if(a==b){h=0;}//ss_sigma
    else if(a==0){h=V[1]*(vconum[b-1]-(*d).at(b-1)*(*d).at(conum))/distr;}//sp_sigma row
    else if(b==0){h=-V[1]*(vconum[a-1]-(*d).at(a-1)*(*d).at(conum))/distr;}//sp_sigma column
  }
  else {h=(V[2]-V[3])*(vconum[a-1]*(*d).at(b-1)+vconum[b-1]*(*d).at(a-1)-2*(*d).at(a-1)*(*d).at(b-1)*(*d).at(conum))/distr;}//pp_sigma and pp_pi off-diagonal

  //V&G routine ends
  return h;
} //Hamder() ends

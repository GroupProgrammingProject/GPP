#include "../include/kpointsfunctions.h"
// Read in kpoints from a file (e.g. for a calculation)
void readinkpoints(char* filename, std::vector<std::pair<std::vector<double>,double> >* kpoints, bool ksymm) {
  std::ifstream infile(filename);
  double k0, k1, k2;
  std::vector<std::vector<double> > kpointstemp;
  std::vector<std::vector<double> > kpointstemp2;
  std::vector<double> kweights;
  std::vector<double> kvec(3);
  std::vector<std::vector<double> >::iterator element;
  int welement;
  typedef std::pair<std::vector<double>,double> pair; 
  while (infile>>k0>>k1>>k2){
	 if (ksymm == 1) {
		kvec.at(0) = std::abs(k0);
		kvec.at(1) = std::abs(k1);
		kvec.at(2) = std::abs(k2);
	 }
	 else {
		kvec.at(0) = k0;
		kvec.at(1) = k1;
		kvec.at(2) = k2;
	 }
	 kpointstemp.push_back(kvec);
  }
  // Have now fully formed kpointstemp object, need to calculate weighting
  int ktot = kpointstemp.size();
  double wfract = 1.0/((double)ktot);
  // Cycle through kvectors checking which are the same by symmetry
  for (int i=0;i<ktot;i++) {
	 // If kpointstemp.at(i) not already in kpointstemp2, sum over j and count instances of it
	 element = std::find(kpointstemp2.begin(), kpointstemp2.end(), kpointstemp.at(i));        // Element in kpointstemp2
	 if(element == kpointstemp2.end()) {                                                      // Check element exists
		for (int j=0;j<ktot;j++) {
		  if ( kpointstemp.at(i)==kpointstemp.at(j) ) {              // Indicates that vectors are equal
			 // If kpointstemp.at(i) not already in kpointstemp2, add it
			 element = std::find(kpointstemp2.begin(), kpointstemp2.end(), kpointstemp.at(i));        // Element in kpointstemp2
			 if (element == kpointstemp2.end()) {
				kpointstemp2.push_back(kpointstemp.at(i));
				kweights.push_back(wfract);
			 }
			 // Otherwise just increase the weighting by one instance
			 else {
				welement = std::distance(kpointstemp2.begin(),element);
				kweights.at(welement) = kweights.at(welement) + wfract;
			 }
		  }
		}
	 }
  }
  int ktot2 = kpointstemp2.size();
  // Now form final object to store kpoints
  for (int i=0;i<ktot2;i++) {
	 (*kpoints).push_back(pair(kpointstemp2.at(i),kweights.at(i)));
  }
}

// Generate a user specified grid of kpoints (only for orthorhombic symmetries)
void genkgrid(char* filename, std::vector<double>* lats, int kgrid[3], bool gamma) {
  std::ofstream outfile(filename);
  bool gammas[3];
  for (int i=0;i<3;i++) {
	 gammas[i] = gamma;                           // Centre on gamma point in each dimension according to input parameter  
	 if (kgrid[i] == 1) {gammas[i] = 1;}          // If only a single point needed make sure to centre on gamma
  }
  std::cout << "kpoint grid " << kgrid[0] << " " << kgrid[1] << " " << kgrid[2] << std::endl;
  double kstep[3] = {2*M_PI/((*lats).at(0)*kgrid[0]), 2*M_PI/((*lats).at(1)*kgrid[1]), 2*M_PI/((*lats).at(2)*kgrid[2])}; 
  double kvec[3];
  for(double k0=-M_PI/(*lats).at(0);k0<M_PI/(*lats).at(0)-0.5*kstep[0];k0=k0+kstep[0]){
	 for(double k1=-M_PI/(*lats).at(1);k1<M_PI/(*lats).at(1)-0.5*kstep[1];k1=k1+kstep[1]){
		for(double k2=-M_PI/(*lats).at(2);k2<M_PI/(*lats).at(2)-0.5*kstep[2];k2=k2+kstep[2]){
		  kvec[0] = k0 + gammas[0]*0.5*kstep[0];
		  kvec[1] = k1 + gammas[1]*0.5*kstep[1];
		  kvec[2] = k2 + gammas[2]*0.5*kstep[2];
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

void kforces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXcd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz, std::vector<double>* kvec, std::vector<double>* TBparam)
{ int k,i,j,l,lp,n,m,nearlabel; /* dummy indeces for cycles*/
  std::vector<double> ddnorm(3);
  double sumphinn,sumphi,derivx,derivy,derivz,dx,dy,dz,r,kr,compderx,compdery,compderz;
  double dualeigen1, dualeigen2,realderx,realdery,realderz,imagderx,imagdery,imagderz;
  double sr,dsrx,dsry,dsrz,gh,ca,cb,cc,cd,skr,ckr;
  Eigen::MatrixXd eigvectr(norbs*N,norbs*N);
  Eigen::MatrixXd eigvecti(norbs*N,norbs*N);
  eigvectr=(*c).real();
  eigvecti=(*c).imag();

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
      if(r<rc){sumphi=sumphi+o(r);}
    }
    for(j=0;j<(*nnear).at(i);j++){ /*Cycle spanning the nearest neighbours of i*/
      nearlabel=(*inear)(i,j);
      sumphinn=0;
      for(m=0;m<(*nnear).at(nearlabel);m++){     
		  dx=(*rx)(nearlabel,(*inear)(nearlabel,m)); /*Definition of vector distances*/
		  dy=(*ry)(nearlabel,(*inear)(nearlabel,m));
		  dz=(*rz)(nearlabel,(*inear)(nearlabel,m));
		  r=(*modr)(nearlabel,(*inear)(nearlabel,m)); /*Modulus of distance*/
		  if(r<rc){sumphinn=sumphinn+o(r);}
      }
      dx=(*rx)(i,nearlabel); /*Definition of vector distances*/
      dy=(*ry)(i,nearlabel);
      dz=(*rz)(i,nearlabel);
      
      r=(*modr)(i,nearlabel); /*Modulus of distance*/
      if(r<rc){
		  ddnorm.at(0)=dx/r; /*Normalised Cartesian vector distances*/
		  ddnorm.at(1)=dy/r;
		  ddnorm.at(2)=dz/r;
		  
		  // Terms which need not be recalculated every orbital iteration
		  kr = dx*(*kvec).at(0) + dy*(*kvec).at(1) + dz*(*kvec).at(2);
		  ckr=cos(kr);
		  skr=sin(kr);
		  sr=s(r);
		  dsrx=ds(r,dx); //derivative of scaling function wrt x
		  dsry=ds(r,dy);
		  dsrz=ds(r,dz);
		  
		  for(l=0;l<norbs;l++){ /*Cycle spanning the first orbital type*/
			 for(lp=0;lp<norbs;lp++){ /*Cycle spanning the second orbital type*/
				gh=Gethijab(i,nearlabel,l,lp,&ddnorm,TBparam); //Hamiltonian matrix element (unscaled)
				derivx=dsrx*gh+sr*kHamder(i,nearlabel,l,lp,&ddnorm,r,0,TBparam);
				derivy=dsry*gh+sr*kHamder(i,nearlabel,l,lp,&ddnorm,r,1,TBparam);
				derivz=dsrz*gh+sr*kHamder(i,nearlabel,l,lp,&ddnorm,r,2,TBparam);

				compderx=(*kvec).at(0)*sr*gh; //term which arises due to complex Hamiltonian
				compdery=(*kvec).at(1)*sr*gh;
				compderz=(*kvec).at(2)*sr*gh;
				
				//split Hamiltonian into real and imaginary matrices: H=R+i*I
				realderx=ckr*derivx-skr*compderx; //dR/dx
				realdery=ckr*derivy-skr*compdery;
				realderz=ckr*derivz-skr*compderz;
				imagderx=skr*derivx+ckr*compderx; //dI/dx
				imagdery=skr*derivy+ckr*compdery;
				imagderz=skr*derivz+ckr*compderz;
				
				for(n=0;n<norbs*N;n++){ /*Cycle spanning the level of the eigenvector*/
				  ca=eigvectr(n,l+i*norbs);
				  cb=eigvecti(n,l+i*norbs);
				  cc=eigvectr(n,lp+nearlabel*norbs);
				  cd=eigvecti(n,lp+nearlabel*norbs);
				  dualeigen1=ca*cc+cb*cd; //c1=ca+i*cb, c2=cc+i*cd
				  dualeigen2=ca*cd-cb*cc;

				  (*fx).at(i)=(*fx).at(i)-2*(realderx*dualeigen1-imagderx*dualeigen2);
				  (*fy).at(i)=(*fy).at(i)-2*(realdery*dualeigen1-imagdery*dualeigen2);
				  (*fz).at(i)=(*fz).at(i)-2*(realderz*dualeigen1-imagderz*dualeigen2);
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
double kHamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum,std::vector<double>* TBparam){
  double h,V[4];//h,Es,Ep and V[4] is only used locally in kHamder()
	V[0]=TBparam->at(2);
	V[1]=TBparam->at(3);
	V[2]=TBparam->at(4);
	V[3]=TBparam->at(5);
	//	CC interaction 0=ss_sigma, 1=sp_sigma, 2=pp_sigma, 3=pp_pi
	// V[0]=-5;V[1]=4.7;V[2]=5.5;V[3]=-1.55;

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

// Finds average energy from multiple kpoints
double avekenergy(int N,int norbs,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr,std::vector<std::pair<std::vector<double>,double> >* kpoints,std::vector<double>* TBparam) {
  int ktot = (*kpoints).size();
  std::vector<double> kvec(3);
  double ebstemp=0;
  double kweight;
  Eigen::MatrixXcd eigvects(4*N,4*N);
  std::vector<double> eigenvalaar (4*N);
  bool v = 0;
  // Get energies over kpoints and average
  for (int i=0;i<ktot;i++) {
	 kvec = (*kpoints).at(i).first;
	 kweight = (*kpoints).at(i).second;
	 std::cout << "k = " << kvec.at(0) << "\t" << kvec.at(1) << "\t" << kvec.at(2) << std::endl;
	 std::cout << "weight = " << kweight << std::endl;
	 ebstemp=ebstemp+kweight*band_Hamiltonian(N,norbs,TBparam,modr,rx,ry,rz,&eigvects,&eigenvalaar,&kvec,v);
	 std::cout << "ebs = " << (1.0/((double)N))*band_Hamiltonian(N,norbs,TBparam,modr,rx,ry,rz,&eigvects,&eigenvalaar,&kvec,v) << std::endl;
  }
  return ebstemp;
}

// Finds average force from multiple kpoints
double avekforces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz, std::vector<std::pair<std::vector<double>,double> >* kpoints,std::vector<double>* TBparam) {
  int ktot = (*kpoints).size();
  std::vector<double> kvec(3);
  double ebstemp=0;
  double kweight;
  Eigen::MatrixXcd eigvects(4*N,4*N);
  std::vector<double> eigenvalaar (4*N);
  std::vector<double> fxtemp(N);
  std::vector<double> fytemp(N);
  std::vector<double> fztemp(N);
  bool v = 0;
  // zero all forces
  for (int l=0;l<N;l++) {
	 (*fx).at(l) = 0;
	 (*fy).at(l) = 0;
	 (*fz).at(l) = 0;
  }

  // Get forces over kpoints and average
  for (int i=0;i<ktot;i++) {
	 kvec = (*kpoints).at(i).first;
	 kweight = (*kpoints).at(i).second;
	 ebstemp=ebstemp+kweight*band_Hamiltonian(N,norbs,TBparam,modr,rx,ry,rz,&eigvects,&eigenvalaar,&kvec,v);
	 kforces(N,4,rc, rx, ry, rz, modr, &eigvects, nnear, inear, &fxtemp, &fytemp, &fztemp, &kvec,TBparam);
	 // Add forces from different kpoints
	 for (int l=0;l<N;l++) {
		(*fx).at(l) = (*fx).at(l) + (1/(double)ktot)*fxtemp.at(l);
		(*fy).at(l) = (*fy).at(l) + (1/(double)ktot)*fytemp.at(l);
		(*fz).at(l) = (*fz).at(l) + (1/(double)ktot)*fztemp.at(l);
	 }
  }
  std::cout << "aveforce ebs = " << ebstemp << std::endl;
  return ebstemp;
}

double kverlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, double &ebs, std::vector<double>* lats, bool pbc,double T, double nu,bool ander, std::vector<std::pair<std::vector<double>,double> >* kpoints,std::vector<double>* TBparam)
{ double boltz=1./11603,svxm=0.0,svym=0.0,svzm=0.0,kin,Tf,sigma=sqrt(boltz*T*m),vxm=0.0,vym=0,vzm=0,rang;
  int N=(*x).size(),nearlabel;
  bool renn=0,v=0;
  std::vector<double> fx(N),fy(N),fz(N),fxn(N),fyn(N),fzn(N);
  Ran ran(time(0)+clock());
  ebs=avekforces(N,norbs,rc,rx,ry,rz,modr,nnear,inear,&fx,&fy,&fz,kpoints,TBparam); //calculate the forces
  //forces(N,norbs,rc,rx,ry,rz,modr,c,nnear,inear,&fx,&fy,&fz,TBparam); //calculate the forces
  for(int i=0; i<N; i++){
    (*x).at(i)=(*x).at(i)+(*vx).at(i)*dt+0.5*fx.at(i)*dt*dt/m;
    (*y).at(i)=(*y).at(i)+(*vy).at(i)*dt+0.5*fy.at(i)*dt*dt/m;
    (*z).at(i)=(*z).at(i)+(*vz).at(i)*dt+0.5*fz.at(i)*dt*dt/m;
  }
  if(pbc==1){pbcshift(x,y,z,lats);} //Put all atoms in unit cell if using PBCs
  renn=RecalculateNearestNeighbours(refx,refy,refz,x,y,z,rc,rv);  
  GetDistances(modr,rx,ry,rz,x,y,z,lats,rv,pbc);
  if(renn==1){
    GetDistances(modr,rx,ry,rz,x,y,z,lats,rv,pbc);
    NearestNeighbours(inear,nnear,modr,rv);
  }
  else{
    for(int i=0;i<N;i++){
      for(int j=0;j<(*nnear).at(i);j++){
	nearlabel=(*inear)(i,j);
	(*rx)(i,nearlabel)=(*x).at(i)-(*x).at(nearlabel);	
	(*ry)(i,nearlabel)=(*y).at(i)-(*y).at(nearlabel);
	(*rz)(i,nearlabel)=(*z).at(i)-(*z).at(nearlabel);
	if(pbc==1){
	  if((*rx)(i,nearlabel)>(*lats).at(0)/2){
	    (*rx)(i,nearlabel)=(*rx)(i,nearlabel)-(*lats).at(0);
	  }
	  else if((*rx)(i,nearlabel)<-(*lats).at(0)/2){
	    (*rx)(i,nearlabel)=(*rx)(i,nearlabel)+(*lats).at(0);
	  }
	  if((*ry)(i,nearlabel)>(*lats).at(1)/2){
	    (*ry)(i,nearlabel)=(*ry)(i,nearlabel)-(*lats).at(1);
	  }
	  else if((*ry)(i,nearlabel)<-(*lats).at(1)/2){
	    (*ry)(i,nearlabel)=(*ry)(i,nearlabel)+(*lats).at(1);
	  }
	  if((*rz)(i,nearlabel)>(*lats).at(2)/2){
	    (*rz)(i,nearlabel)=(*rz)(i,nearlabel)-(*lats).at(2);
	  }
	  else if((*rz)(i,nearlabel)<-(*lats).at(2)/2){
	    (*rz)(i,nearlabel)=(*rz)(i,nearlabel)+(*lats).at(2);
	  }
	}
	(*modr)(i,nearlabel)=sqrt((*rx)(i,nearlabel)*(*rx)(i,nearlabel)+(*ry)(i,nearlabel)*(*ry)(i,nearlabel)+(*rz)(i,nearlabel)*(*rz)(i,nearlabel));
      }
    }
  }
  ebs=avekforces(N,norbs,rc,rx,ry,rz,modr,nnear,inear,&fxn,&fyn,&fzn,kpoints,TBparam); //recalculate the forces
  //ebs=Hamiltonian(N,norbs,TBparam,modr,rx,ry,rz,c,v);
  //forces(N,norbs,rc,rx,ry,rz,modr,c,nnear,inear,&fxn,&fyn,&fzn,TBparam);//recalculate forces
  for(int i=0; i<N; i++)//calculate new velocities
    {
      (*vx).at(i)=(*vx).at(i)+dt*(fx.at(i)+fxn.at(i))/(2*m);
      (*vy).at(i)=(*vy).at(i)+dt*(fy.at(i)+fyn.at(i))/(2*m);
      (*vz).at(i)=(*vz).at(i)+dt*(fz.at(i)+fzn.at(i))/(2*m);
      rang=ran.doub();
      if((rang<nu*dt) && (ander==1)){ //implement the Andersen thermostat for canonical ensemble
	(*vx).at(i)=Gauss(0,sigma)/m; //generate random numbers from Gaussian distribution
	(*vy).at(i)=Gauss(0,sigma)/m;
	(*vz).at(i)=Gauss(0,sigma)/m;
	vxm=vxm+(*vx).at(i);
	vym=vym+(*vy).at(i);
	vzm=vzm+(*vz).at(i);
      }
    }
  for(int i=0; i<N; i++)//mean square velocities
    {
      (*vx).at(i)=(*vx).at(i)-vxm/N;
      (*vy).at(i)=(*vy).at(i)-vym/N;
      (*vz).at(i)=(*vz).at(i)-vzm/N;
      svxm=svxm+(*vx).at(i)*(*vx).at(i);
      svym=svym+(*vy).at(i)*(*vy).at(i);
      svzm=svzm+(*vz).at(i)*(*vz).at(i);
    }
  kin=0.5*m*(svxm+svym+svzm); //kinetic energy
  Tf=2*kin/(3*boltz*(N-1)); //final temperature  
  return Tf;	
}



#include "../include/MolDyn.h"

/*Inputs, in order: #orbitals (norbs); cut-off radius (rc); Verlet radius (rv); atomic mass (m); timestep [ps] (dt); position arrays (x,y,z); reference position arrays (refx,refy,refz); velocity arrays (vx,vy,vz); matrix of eigenvectors, with vectors as rows (c); nearest neighbour lists (vector nnear and matrix inear); matrices of vector distances (rx,ry,rz); matrix of moduli of vector distances (modr); band structure energy (ebs); supercell dimensions (lats); periodic boundary conditions switch (pbc); temperature (T); frequency of velocity rescalings for thermostat (nu); thermostat switch (ander); tight binding parameters (TBparam).*/
double verlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, Eigen::MatrixXd* c,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr,double &ebs, std::vector<double>* lats, bool pbc,double T, double nu,bool ander,std::vector<double>* TBparam)
{ double boltz=1./11603,svxm=0.0,svym=0.0,svzm=0.0,kin,Tf,sigma=sqrt(boltz*T*m),vxm=0.0,vym=0,vzm=0,rang;
  int N=(*x).size(),nearlabel;
  bool renn=0,v=0;
  std::vector<double> fx(N),fy(N),fz(N),fxn(N),fyn(N),fzn(N);
  /*initialisation of seed for RNG:*/
  Ran ran(time(0)+clock());
  /*call to forces function:*/
  forces(N,norbs,rc,rx,ry,rz,modr,c,nnear,inear,&fx,&fy,&fz,TBparam);
  /*integration of atomic positions:*/
  for(int i=0; i<N; i++){
    (*x).at(i)=(*x).at(i)+(*vx).at(i)*dt+0.5*fx.at(i)*dt*dt/m;
    (*y).at(i)=(*y).at(i)+(*vy).at(i)*dt+0.5*fy.at(i)*dt*dt/m;
    (*z).at(i)=(*z).at(i)+(*vz).at(i)*dt+0.5*fz.at(i)*dt*dt/m;
  }
  /*if pbc are active, all atoms are shifted to the supercell:*/
  if(pbc==1){pbcshift(x,y,z,lats);}
  /*call to RecalculateNearestNeighbours:*/
  renn=RecalculateNearestNeighbours(refx,refy,refz,x,y,z,rc,rv);
  /*if renn==1 (one atom moved with respect to it reference position more than the threshold), recalculate NNs:*/
  if(renn==1){
    GetDistances(modr,rx,ry,rz,x,y,z,lats,rv,pbc);
    NearestNeighbours(inear,nnear,modr,rv);
  }
  /*else, perform distance calculations only on NNs:*/
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
  /*call to Hamiltonian for new eigenvectors*/
  ebs=Hamiltonian(N,norbs,TBparam,modr,rx,ry,rz,c,v);
  /*recalculation of forces:*/
  forces(N,norbs,rc,rx,ry,rz,modr,c,nnear,inear,&fxn,&fyn,&fzn,TBparam);
  /*integration of velocities:*/
  for(int i=0; i<N; i++)
    {
      (*vx).at(i)=(*vx).at(i)+dt*(fx.at(i)+fxn.at(i))/(2*m);
      (*vy).at(i)=(*vy).at(i)+dt*(fy.at(i)+fyn.at(i))/(2*m);
      (*vz).at(i)=(*vz).at(i)+dt*(fz.at(i)+fzn.at(i))/(2*m);
      rang=ran.doub();
      /*if thermostat is on, produce new Maxwell-distributed velocities, with frequency nu:*/
      if((rang<nu*dt) && (ander==1)){ 
	(*vx).at(i)=Gauss(0,sigma)/m; 
	(*vy).at(i)=Gauss(0,sigma)/m;
	(*vz).at(i)=Gauss(0,sigma)/m;
	vxm=vxm+(*vx).at(i);
	vym=vym+(*vy).at(i);
	vzm=vzm+(*vz).at(i);
      }
    }
  /*calculation of mean squared velocities:*/
  for(int i=0; i<N; i++)
    {
		if(ander==1){
			(*vx).at(i)=(*vx).at(i)-vxm/N;
			(*vy).at(i)=(*vy).at(i)-vym/N;
			(*vz).at(i)=(*vz).at(i)-vzm/N;
		}
      svxm=svxm+(*vx).at(i)*(*vx).at(i);
      svym=svym+(*vy).at(i)*(*vy).at(i);
      svzm=svzm+(*vz).at(i)*(*vz).at(i);
    }
  /*calculation of kinetic energy and final temperature:*/
  kin=0.5*m*(svxm+svym+svzm); 
  Tf=2*kin/(3*boltz*(N-1)); 
  return Tf;	
}

/*Inputs, in order: #atoms (N); #orbitals (norbs); cut-off radius (rc); matrices of vector distances (rx,ry,rz); matrix of moduli of vector distances (modr); matrix of eigenvectors, with vectors as rows (c);  nearest neighbour lists (vector nnear and matrix inear); forces vectors (fx,fy,fz); tight binding parameters (TBparam).*/
void forces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz,std::vector<double>* TBparam)
{ int k,i,j,l,lp,n,m,nearlabel;
  std::vector<double> ddnorm(3);
  double sumphinn,sumphi,dualeigen,derivx,derivy,derivz,dx,dy,dz,r,dsrx,dsry,dsrz,sr,gh;
  /*initialisation of forces:*/
  for(i=0;i<N;i++){
    (*fx).at(i)=0;
    (*fy).at(i)=0;
    (*fz).at(i)=0;
  }
  /*cycle to compute forces on atom i:*/
  for(i=0;i<N;i++){
    sumphi=0;
    /*calculation of sumphi, i.e. sum of pairwise phi(r_ij) over NNs of i. Needed for repulsive forces.*/
    for(k=0;k<(*nnear).at(i);k++){
      nearlabel=(*inear)(i,k);
      r=(*modr)(i,nearlabel);
      if(r<rc){
	sumphi=sumphi+o(r);
      }
    }
    /*cycle over NNs of i*/
    for(j=0;j<(*nnear).at(i);j++){
      nearlabel=(*inear)(i,j);
      sumphinn=0;
      /*calculation of sumphinn, i.e. sum of pairwise phi(r_jk) over NNs of j. Needed for repulsive forces.*/
      for(m=0;m<(*nnear).at(nearlabel);m++){     
	dx=(*rx)(nearlabel,(*inear)(nearlabel,m));
	dy=(*ry)(nearlabel,(*inear)(nearlabel,m));
	dz=(*rz)(nearlabel,(*inear)(nearlabel,m));
	r=(*modr)(nearlabel,(*inear)(nearlabel,m));
	if(r<rc){
	  sumphinn=sumphinn+o(r);
	}
      }
      /*exctracting relevant distances from the distance matrices:*/
      dx=(*rx)(i,nearlabel);
      dy=(*ry)(i,nearlabel);
      dz=(*rz)(i,nearlabel); 
      r=(*modr)(i,nearlabel);
      /*if |r_ij|<rc, compute forces:*/
      if(r<rc){
	/*definition of normalised distance vectors:*/
	ddnorm.at(0)=dx/r;
	ddnorm.at(1)=dy/r;
	ddnorm.at(2)=dz/r;
	/*call to ds to compute derivatives of scaling function:*/
	dsrx=ds(r,dx);
	dsry=ds(r,dy);
	dsrz=ds(r,dz);
	sr=s(r);
	/*cycle over the first orbital type:*/
	for(l=0;l<norbs;l++){
	  /*cycle over the second orbital type:*/
	  for(lp=0;lp<norbs;lp++){ 
	    /*call to Gethijab to get relevant TB matrix element*/
	    gh=Gethijab(i,nearlabel,l,lp,&ddnorm,TBparam);
	    /*calculation of derivatives of matrix element:*/
	    derivx=ds(r,dx)*gh+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,0,TBparam);
	    derivy=ds(r,dy)*gh+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,1,TBparam);
	    derivz=ds(r,dz)*gh+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,2,TBparam);
	    /*cycle over eigenvector numbers:*/
	    for(n=0;n<norbs*N;n++){
	      /*definition of product of relevant eigenvectors:*/
	      dualeigen=(*c)(n,l+i*norbs)*(*c)(n,lp+nearlabel*norbs);  
	      /*calculation of band structure forces:*/
	      (*fx).at(i)=(*fx).at(i)-2*derivx*dualeigen;
	      (*fy).at(i)=(*fy).at(i)-2*derivy*dualeigen;
	      (*fz).at(i)=(*fz).at(i)-2*derivz*dualeigen;
	    }
	  }
	}
	/*calculation of repulsive forces:*/
	(*fx).at(i)=(*fx).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dx);
	(*fy).at(i)=(*fy).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dy);
	(*fz).at(i)=(*fz).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dz);
      }
    }
  }  
}

/*Inputs, in order: atomic mass (m); velocity arrays (vx,vy,vz); initial temperature (T).*/
void velocity(double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T)
{ int N=(*vx).size();
  double c,ke,boltz=1./11603;//1.38*pow(10,-23);
  double vxtot=0,vytot=0,vztot=0,msvx=0,msvy=0,msvz=0,vxavg,vyavg,vzavg,Tp;
  /*definition of bounds of uniform distribution (to get average and standard deviation equal to Maxwell's distribution):*/
  c=sqrt(3*boltz*T/m);
  /*initialisation of seed for RNG:*/
  Ran ran(time(0));
  /*cycle to compute velocities:*/
  for (int i=0; i<N; i++)
    {  
      (*vx).at(i)=c*(2*ran.doub()-1);
      (*vy).at(i)=c*(2*ran.doub()-1);
      (*vz).at(i)=c*(2*ran.doub()-1);
      vxtot=vxtot+(*vx).at(i);
      vytot=vytot+(*vy).at(i);
      vztot=vztot+(*vz).at(i);
    }
  /*calculation of average velocities:*/
  vxavg=vxtot/N;
  vyavg=vytot/N;
  vzavg=vztot/N;
  /*velocity rescaling to ensure zero net momentum:*/
  for (int i=0; i<N; i++)
    {
      (*vx).at(i)=(*vx).at(i)-vxavg;
      (*vy).at(i)=(*vy).at(i)-vyavg;
      (*vz).at(i)=(*vz).at(i)-vzavg;
      msvx=msvx+(*vx).at(i)*(*vx).at(i);
      msvy=msvy+(*vy).at(i)*(*vy).at(i);
      msvz=msvz+(*vz).at(i)*(*vz).at(i);
    }
  /*calculation of kinetic energy:*/
  ke=0.5*m*(msvx+msvy+msvz);
  /*calculation of temperature*/
  Tp=2*ke/(3*boltz*(N-1));
  /*velocity rescaling to ensure actual prescribed temperature:*/
  for(int i=0; i<N; i++)
    {
      (*vx).at(i)=(*vx).at(i)*sqrt(T/Tp);
      (*vy).at(i)=(*vy).at(i)*sqrt(T/Tp);
      (*vz).at(i)=(*vz).at(i)*sqrt(T/Tp);
    }
}

/*Inputs, in order: atomic indeces (i,j); orbital indeces (a,b); distance vector between i and j (d); modulus of distance between i and j (distr); required derivative component (conum, i.e. 0=x, 1=y, 2=z); Tight Binding parameters (TBparam).*/
double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum,std::vector<double>* TBparam){
  double h,V[4];
  /*retrieving relevant TB parameters:*/
  V[0]=TBparam->at(2);
  V[1]=TBparam->at(3);
  V[2]=TBparam->at(4);
  V[3]=TBparam->at(5);
  /*definition of a vector with only non zero component vconum[conum]. Used to implement a representation of Kronecker delta.*/
  int vconum[3]; 
  vconum[0]=0;
  vconum[1]=0;
  vconum[2]=0;
  vconum[conum]=1;
  /*if statements to reproduce the correct form of the required derivative:*/
  if(i==j){h=0;}
  else if(a*b==0){
    if(a==b){h=0;}//ss_sigma
    else if(a==0){h=V[1]*(vconum[b-1]-(*d).at(b-1)*(*d).at(conum))/distr;}//sp_sigma row
    else if(b==0){h=-V[1]*(vconum[a-1]-(*d).at(a-1)*(*d).at(conum))/distr;}//sp_sigma column
  }
  else {h=(V[2]-V[3])*(vconum[a-1]*(*d).at(b-1)+vconum[b-1]*(*d).at(a-1)-2*(*d).at(a-1)*(*d).at(b-1)*(*d).at(conum))/distr;}
  /*returning the derivative of the required element:*/
  return h;
}

/*Inputs, in order: #orbitals (norbs); cut-off radius (rc); Verlet radius (rv); atomic mass (m); timestep [ps] (dt); position arrays (x,y,z); reference position arrays (refx,refy,refz); matrix of eigenvectors, with vectors as rows (eigvects); nearest neighbour lists (vector nnear and matrix inear); matrices of vector distances (rx,ry,rz); matrix of moduli of vector distances (modr); supercell dimensions (lats); periodic boundary conditions switch (pbc); temperature (T); frequency of velocity rescalings for thermostat (nu); integration step for steepest descent (h); output switch (verb); interval of steps printed to file (nprint); tight binding parameters (TBparam); steepest descent tolerance on max of moduli of forces (tol); maximum number of steepest descent steps (maxsteep); k-points switch (kpts); kpoint grid (kpoints).*/
int GeomOpt(int norbs,double rc,double rv,double m,double dt,int nmd,std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* refposx, std::vector<double>* refposy, std::vector<double>* refposz, Eigen::MatrixXd* eigvects,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, std::vector<double>* lats, bool pbc,double T,double nu,double h,bool verb, int nprint,std::vector<double>* TBparam, double tol,int maxsteep,bool kpts,std::vector<std::pair<std::vector<double>,double> >* kpoints){
  int n=(*posx).size();
  if (kpts==1) {std::cout << "ktot = " << (*kpoints).size() << std::endl;}
  bool v=0, renn=0, ander=1;
  int i,j,nearlabel,nsteep;
  double Tin=T,Tf,tmd,kb=1./11603;
  std::vector<double> vx(n), vy(n), vz(n), fx(n), fy(n), fz(n);
  double ebs,erep,etot,ekin,fmax;
  /*calculation of initial velocities:*/
  velocity(m,&vx,&vy,&vz,T);
  /*opening of output files:*/
  FILE *file=fopen("movie_relax.txt","w");
  FILE *en=fopen("energy_relax.txt","w");
  FILE *file2=fopen("forces_relax.txt","w");
  /*if verb=1, print initial positions and energy to file*/
  if(verb==1){
    fprintf(file,"%d\nIteration %d\n",n,i);
    for(i=0;i<n;i++){
      fprintf(file,"6  %f %f %f\n",(*posx).at(i),(*posy).at(i),(*posz).at(i));
    }
    /*calculation of initial band structure energy and eigenvectors:*/
    if (kpts==1) {ebs=avekenergy(n,norbs,rx,ry,rz,modr,kpoints,TBparam);}
    else {ebs=Hamiltonian(n,norbs,TBparam,modr,rx,ry,rz,eigvects,v);}
    /*calculation of repulsive, kinetic and total energy:*/
    erep=Erep(modr);	
    ekin=3*(n-1)*kb*T/2;
    etot=ebs+erep+ekin;
    fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",tmd,Tf,ekin,ebs,erep,etot);
  }
  std::cout << "Starting simulated annealing..." << std::endl;
  /*simulated annealing cycle:*/
  for(i=1;i<nmd+1;i++){
    if(i%(nmd/10)==0){std::cout << i*100/nmd << "% completed" << std::endl;}
    /*linear decrease of temperature (0.99 factor implies that a longer time is spent at the lowest temperature):*/
    T=T-Tin/(nmd*0.99);
    /*if T gets below zero, redefine it to be T=1e-6:*/
    if(T<=0){
      T=1e-6;
    }
    /*call to verlet algorithm:*/
    if (kpts==1) {Tf=kverlet(norbs,rc,rv,m,dt,posx,posy,posz,refposx,refposy,refposz,&vx,&vy,&vz,nnear,inear,rx,ry,rz,modr,ebs,lats,pbc,T,nu,ander,kpoints,TBparam);}
    else {Tf=verlet(norbs,rc,rv,m,dt,posx,posy,posz,refposx,refposy,refposz,&vx,&vy,&vz,eigvects,nnear,inear,rx,ry,rz,modr,ebs,lats,pbc,T,nu,ander,TBparam);}
    /*printing positions and energies to file (every nprint steps and if verb=1):*/
    if(i%nprint==0 && verb==1){
      erep=Erep(modr);
      etot=ebs+erep+ekin;
      tmd=i*dt;
      fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",tmd,Tf,ekin,ebs,erep,etot);
      fprintf(file,"%d\nIteration %d\n",n,i);
      for(j=0;j<n;j++){
	fprintf(file,"6  %f %f %f\n",(*posx).at(j),(*posy).at(j),(*posz).at(j));
      }
    } 
  }
  /*calculation of new forces*/
  if (kpts==1) {ebs=avekforces(n,norbs,rc,rx,ry,rz,modr,nnear,inear,&fx,&fy,&fz,kpoints,TBparam);}
  else {
    ebs=Hamiltonian(n,norbs,TBparam,modr,rx,ry,rz,eigvects,v);
    forces(n,norbs,rc,rx,ry,rz,modr,eigvects,nnear,inear,&fx,&fy,&fz,TBparam);
  }
  /*calculation of maximum of force moduli:*/
  fmax=0;
  for(j=0;j<n;j++){
    if(fabs(fx.at(j))>fmax){fmax=fabs(fx.at(j));}
    if(fabs(fy.at(j))>fmax){fmax=fabs(fy.at(j));}
    if(fabs(fz.at(j))>fmax){fmax=fabs(fz.at(j));}
  }
  nsteep=0;
  std::cout << "Starting steepest descent..." << endl;
  /*steepest descent cycle:*/
  while(fmax>tol && nsteep<maxsteep){
    fmax=0;
    renn=RecalculateNearestNeighbours(refposx,refposy,refposz,posx,posy,posz,rc,rv); 
    if(renn==1){
      GetDistances(modr,rx,ry,rz,posx,posy,posz,lats,rv,pbc);
      NearestNeighbours(inear,nnear,modr,rv);
    }
    else{
      for(i=0;i<n;i++){
	for(j=0;j<(*nnear).at(i);j++){
	  nearlabel=(*inear)(i,j);
	  (*rx)(i,nearlabel)=(*posx).at(i)-(*posx).at(nearlabel);
	  (*ry)(i,nearlabel)=(*posy).at(i)-(*posy).at(nearlabel);
	  (*rz)(i,nearlabel)=(*posz).at(i)-(*posz).at(nearlabel);
	  if(pbc==1){
	    if((*rx)(i,nearlabel)>(*lats).at(0)/2){
	      (*rx)(i,nearlabel)=(*rx)(i,nearlabel)-(*lats).at(0);
	    }
	    if((*rx)(i,nearlabel)<-(*lats).at(0)/2){
	      (*rx)(i,nearlabel)=(*rx)(i,nearlabel)+(*lats).at(0);
	    }
	    if((*ry)(i,nearlabel)>(*lats).at(1)/2){
	      (*ry)(i,nearlabel)=(*ry)(i,nearlabel)-(*lats).at(1);
	    }
	    if((*ry)(i,nearlabel)<-(*lats).at(1)/2){
	      (*ry)(i,nearlabel)=(*ry)(i,nearlabel)+(*lats).at(1);
	    }
	    if((*rz)(i,nearlabel)>(*lats).at(2)/2){
	      (*rz)(i,nearlabel)=(*rz)(i,nearlabel)-(*lats).at(2);
	    }
	    if((*rz)(i,nearlabel)<-(*lats).at(2)/2){
	    (*rz)(i,nearlabel)=(*rz)(i,nearlabel)+(*lats).at(2);
	    }
	  }
	  (*modr)(i,nearlabel)=sqrt((*rx)(i,nearlabel)*(*rx)(i,nearlabel)+(*ry)(i,nearlabel)*(*ry)(i,nearlabel)+(*rz)(i,nearlabel)*(*rz)(i,nearlabel));
	}
      }
    }
    if (kpts==1) {ebs=avekforces(n,norbs,rc,rx,ry,rz,modr,nnear,inear,&fx,&fy,&fz,kpoints,TBparam);}
    else {
      ebs=Hamiltonian(n,norbs,TBparam,modr,rx,ry,rz,eigvects,v);
      forces(n,norbs,rc,rx,ry,rz,modr,eigvects,nnear,inear,&fx,&fy,&fz,TBparam);
    }
    /*moving atoms positions along forces, with magnitude h*/
    for(j=0;j<n;j++){
      if(fabs(fx.at(j))>fmax){fmax=fabs(fx.at(j));}
      if(fabs(fy.at(j))>fmax){fmax=fabs(fy.at(j));}
      if(fabs(fz.at(j))>fmax){fmax=fabs(fz.at(j));}
      (*posx).at(j)=(*posx).at(j)+h*fx.at(j);
      (*posy).at(j)=(*posy).at(j)+h*fy.at(j);
      (*posz).at(j)=(*posz).at(j)+h*fz.at(j);
    }
    nsteep++;
  }
  GetDistances(modr,rx,ry,rz,posx,posy,posz,lats,rv,pbc);
  /*calculation of final energies:*/
  erep=Erep(modr);
  if (kpts==1) {ebs=avekforces(n,norbs,rc,rx,ry,rz,modr,nnear,inear,&fx,&fy,&fz,kpoints,TBparam);}
  else {ebs=Hamiltonian(n,norbs,TBparam,modr,rx,ry,rz,eigvects,v);}
  etot=ebs+erep;
  /*printing energies and final positions to file if verb==1:*/
  if(verb==1){	
    fprintf(file,"%d\nIteration %d\n",n,i);
    for(j=0;j<n;j++){
      fprintf(file,"6  %f %f %f\n",(*posx).at(j),(*posy).at(j),(*posz).at(j));
      fprintf(file2,"6 %.10f %.10f %.10f \n",fx.at(j),fy.at(j),fz.at(j));
    }
    fprintf(en,"%f\t%f\t%f\t%f\t%f\t%f\n",tmd+1,Tf,ekin,ebs,erep,etot);
  }
  /*output to screen:*/
  std::cout << "Relaxation done." << std::endl;
  std::cout << "Maximum force = " << fmax << std::endl;
  std::cout << "Ebs = " << ebs << std::endl;
  std::cout << "Erep = " << erep << std::endl;
  std::cout << "Etot = " << etot << std::endl;
  
  return 0;
}

#include "../include/MolDyn.h"

double verlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, Eigen::MatrixXd* c,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr,double &ebs, std::vector<double>* lats, bool pbc,double T, double nu,bool ander)
{ double boltz=1./11603,svxm=0.0,svym=0.0,svzm=0.0,kin,Tf,sigma=sqrt(boltz*T*m),vxm=0.0,vym=0,vzm=0,rang;
  int N=(*x).size();
  bool renn=0,v=0;
  std::vector<double> fx(N),fy(N),fz(N),fxn(N),fyn(N),fzn(N);
  Ran ran(time(0)+clock());
  forces(N,norbs,rc,rx,ry,rz,modr,c,nnear,inear,&fx,&fy,&fz); //calculate the forces
  for(int i=0; i<N; i++)
    {
      (*x).at(i)=(*x).at(i)+(*vx).at(i)*dt+0.5*fx.at(i)*dt*dt/m;
      (*y).at(i)=(*y).at(i)+(*vy).at(i)*dt+0.5*fy.at(i)*dt*dt/m;
      (*z).at(i)=(*z).at(i)+(*vz).at(i)*dt+0.5*fz.at(i)*dt*dt/m;
    }
  renn=RecalculateNearestNeighbours(refx,refy,refz,x,y,z,rc,rv);  
  if(renn==1){
    NearestNeighbours(inear,nnear,modr,rv);
  }
  GetDistances(modr,rx,ry,rz,x,y,z,lats,rv,pbc);
  ebs=Hamiltonian(N,modr,rx,ry,rz,c,v);
  forces(N,norbs,rc,rx,ry,rz,modr,c,nnear,inear,&fxn,&fyn,&fzn);//recalculate forces
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
		}
		vxm=vxm+(*vx).at(i);
		vym=vym+(*vy).at(i);
		vzm=vzm+(*vz).at(i);
  }
  vxm=vxm/N;
  vym=vym/N;
  vzm=vzm/N;
  for(int i=0; i<N; i++)//mean square velocities
    {
		(*vx).at(i)=(*vx).at(i)-vxm;
		(*vy).at(i)=(*vy).at(i)-vym;
		(*vz).at(i)=(*vz).at(i)-vzm;
      svxm=svxm+(*vx).at(i)*(*vx).at(i);
      svym=svym+(*vy).at(i)*(*vy).at(i);
      svzm=svzm+(*vz).at(i)*(*vz).at(i);
    }
  kin=0.5*m*(svxm+svym+svzm); //kinetic energy
  Tf=2*kin/(3*boltz*(N-1)); //final temperature 
  return Tf;	
}

/*INPUTS:N=numb. atoms; x,y,z=atom positions (arrays); c=eigenvectors (matrix with each column as the n-th eigenvector); rc=cut-off radius;
nnear=number of nearest neighbours (nn) to i-th atom (array); inear=label of j-th nn to i-th atom (matrix); fx,fy,fz forces on each atom (arrays);
maxnn=max number of nn */
void forces(int N,int norbs,double rc,Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr, Eigen::MatrixXd* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz)
{ int k,i,j,l,lp,n,m,nearlabel; /* dummy indeces for cycles*/
  std::vector<double> ddnorm(3);
  double sumphinn,sumphi,dualeigen,derivx,derivy,derivz,dx,dy,dz,r;

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
	    derivx=ds(r,dx)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,0);
	    derivy=ds(r,dy)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,1);
	    derivz=ds(r,dz)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,2);
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

void velocity(double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T)
{  //calculate velocities
  int N=(*vx).size();
  double c,ke,boltz=1./11603;//1.38*pow(10,-23);
  double vxtot=0,vytot=0,vztot=0,msvx=0,msvy=0,msvz=0,vxavg,vyavg,vzavg,Tp;
  c=sqrt(3*boltz*T/m);
  Ran ran(time(0)+clock()); //CHANGE RNG!!!!!!!!
  for (int i=0; i<N; i++)
    {  
      (*vx).at(i)=c*(2*ran.doub()-1);
      (*vy).at(i)=c*(2*ran.doub()-1);
      (*vz).at(i)=c*(2*ran.doub()-1);
      vxtot=vxtot+(*vx).at(i);
      vytot=vytot+(*vy).at(i);
      vztot=vztot+(*vz).at(i);
    }
  vxavg=vxtot/N;
  vyavg=vytot/N;
  vzavg=vztot/N;
  
  for (int i=0; i<N; i++)
    {
      (*vx).at(i)=(*vx).at(i)-vxavg;
      (*vy).at(i)=(*vy).at(i)-vyavg;
      (*vz).at(i)=(*vz).at(i)-vzavg;
      msvx=msvx+(*vx).at(i)*(*vx).at(i);
      msvy=msvy+(*vy).at(i)*(*vy).at(i);
      msvz=msvz+(*vz).at(i)*(*vz).at(i);
    }
  ke=0.5*m*(msvx+msvy+msvz); //must be adapted for different masses
  Tp=2*ke/(3*boltz*(N-1));
  for(int i=0; i<N; i++)
    {
      (*vx).at(i)=(*vx).at(i)*sqrt(T/Tp);
      (*vy).at(i)=(*vy).at(i)*sqrt(T/Tp);
      (*vz).at(i)=(*vz).at(i)*sqrt(T/Tp);
    }
}

//Hamder() returns value of Hamilatonian matrix element differentiated wrt x,y or z.
double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum){
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

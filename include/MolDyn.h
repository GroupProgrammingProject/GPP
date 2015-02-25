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
#include "functions.h"
#include "Gethijab.h"
#include "geometryinfo.h"
#include "hamiltonian.h"
#include <vector>

//Vector syntax: std::vector<double>. Assignments: xi=(*v).at(i). Inputs as pointer: std::vector<double>*

/*Change vector syntax as: std::vector<double>. Assignments xi=(*v).at(i). Inputs as pointer: std::vector<double>*, to call in functs ad &*/

double verlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>*vy, std::vector<double>* vz, Eigen::MatrixXd* c,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr);/*Inputs, in order: #orbitals; timestep; xyz arrays; velocity arrays; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; atom vector distances; modulus of vector distances.*/

void forces(int N, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz);

void velocity(double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T);/*Inputs, in order: #atoms; mass; vectors of velocities; temperature.*/

double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum);/*Inputs, in order: first atom index; second atom index; first orbital index; second orbital index; vector of direction cosines; modulus of distance between i and j; switch for component.*/ 

double verlet(int norbs,double rc,double rv,double m,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* refx, std::vector<double>* refy, std::vector<double>* refz,std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, Eigen::MatrixXd* c,std::vector<int>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, Eigen::MatrixXd* modr)
{ double boltz=1./11603,svxm=0.0,svym=0.0,svzm=0.0,kin,Tf;
  int N=(*x).size();
  bool renn=0;
  std::vector<double> fx(N),fy(N),fz(N),fxn(N),fyn(N),fzn(N),vxm(N),vym(N),vzm(N);
  forces(N,x,y,z,&fx,&fy,&fz); //calculate the forces
  for(int i=0; i<N; i++)
    {
      (*x).at(i)=(*x).at(i)+(*vx).at(i)*dt+0.5*fx.at(i)*dt*dt/m;
      (*y).at(i)=(*y).at(i)+(*vy).at(i)*dt+0.5*fy.at(i)*dt*dt/m;
      (*z).at(i)=(*z).at(i)+(*vz).at(i)*dt+0.5*fz.at(i)*dt*dt/m;
    }
  renn=RecalculateNearestNeighbours(refx,refy,refz,x,y,z,rc,rv);  if(renn==1){
    NearestNeighbours(inear,nnear,modr,rv);
  }
  GetAllDistances(modr,rx,ry,rz,x,y,z);
  forces(N,x,y,z,&fxn,&fyn,&fzn);//recalculate forces
  for(int i=0; i<N; i++)//calculate new velocities
    {
      (*vx).at(i)=(*vx).at(i)+dt*(fx.at(i)+fxn.at(i))/(2*m);
      (*vy).at(i)=(*vy).at(i)+dt*(fy.at(i)+fyn.at(i))/(2*m);
      (*vz).at(i)=(*vz).at(i)+dt*(fz.at(i)+fzn.at(i))/(2*m);
    }
  for(int i=0; i<N; i++)//mean square velocities
    {
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
void forces(int n, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz)
{ int k,i,j,l,lp,m,nearlabel; /* dummy indeces for cycles*/
  std::vector<double> ddnorm(3);
  double r,h=0.0001;
  Eigen::MatrixXd modrl(n,n);
  Eigen::MatrixXd dlrx(n,n);
  Eigen::MatrixXd dlry(n,n);
  Eigen::MatrixXd dlrz(n,n);
  Eigen::MatrixXd modrr(n,n);
  Eigen::MatrixXd drrx(n,n);
  Eigen::MatrixXd drry(n,n);
  Eigen::MatrixXd drrz(n,n);
  Eigen::MatrixXd eigvects(4*n,4*n);
  std::vector<double> drposx(n),drposy(n),drposz(n),dlposx(n),dlposy(n),dlposz(n);
  bool v=0;

  for(i=0;i<n;i++){ /*initialisation of forces*/
    (*fx).at(i)=0;
    (*fy).at(i)=0;
    (*fz).at(i)=0;
  }
  
  for(i=0;i<n;i++){ /*Cycle to compute band structure forces on atom i*/
    for(j=0;j<n;j++){
      drposx.at(j)=(*posx).at(j);
      dlposx.at(j)=(*posx).at(j);
      drposy.at(j)=(*posy).at(j);
      dlposy.at(j)=(*posy).at(j);
      drposz.at(j)=(*posz).at(j);
      dlposz.at(j)=(*posz).at(j);
    }
    drposx.at(i)=(*posx).at(i)+h;
    dlposx.at(i)=(*posx).at(i)-h;
    GetAllDistances(&modrr,&drrx,&drry,&drrz,&drposx,posy,posz);
    GetAllDistances(&modrl,&dlrx,&dlry,&dlrz,&dlposx,posy,posz);
    (*fx).at(i)=(*fx).at(i)-(Hamiltonian(n,&modrr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&modrl,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
    (*fx).at(i)=(*fx).at(i)-(Erep(&modrr)-Erep(&modrl))/(2*h);
    
    drposy.at(i)=(*posy).at(i)+h;
    dlposy.at(i)=(*posy).at(i)-h;
    GetAllDistances(&modrr,&drrx,&drry,&drrz,posx,&drposy,posz);
    GetAllDistances(&modrl,&dlrx,&dlry,&dlrz,posx,&dlposy,posz);
    (*fy).at(i)=(*fy).at(i)-(Hamiltonian(n,&modrr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&modrl,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
    (*fy).at(i)=(*fy).at(i)-(Erep(&modrr)-Erep(&modrl))/(2*h);
    
    drposz.at(i)=(*posz).at(i)+h;
    dlposz.at(i)=(*posz).at(i)-h;
    GetAllDistances(&modrr,&drrx,&drry,&drrz,posx,posy,&drposz);
    GetAllDistances(&modrl,&dlrx,&dlry,&dlrz,posx,posy,&dlposz);
    (*fz).at(i)=(*fz).at(i)-(Hamiltonian(n,&modrr,&drrx,&drry,&drrz,&eigvects,v)-Hamiltonian(n,&modrl,&dlrx,&dlry,&dlrz,&eigvects,v))/(2*h);
    (*fz).at(i)=(*fz).at(i)-(Erep(&modrr)-Erep(&modrl))/(2*h);
  }
}


void velocity(double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T)
{  //calculate velocities
  int N=(*vx).size();
  double c,ke,boltz=1./11603;//1.38*pow(10,-23);
  double vxtot=0,vytot=0,vztot=0,msvx=0,msvy=0,msvz=0,vxavg,vyavg,vzavg,Tp;
  c=sqrt(3*boltz*T/m);
  srand(time(0)); //CHANGE RNG!!!!!!!!
  for (int i=0; i<N; i++)
    {  
      (*vx).at(i)=c*(2*rand()-1);
      (*vy).at(i)=c*(2*rand()-1);
      (*vz).at(i)=c*(2*rand()-1);
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
  else if((a==b) && (a*b!=0)){h=2*(V[2]-V[3])*(*d).at(a-1)*(vconum[a-1]-(*d).at(a-1)*(*d).at(conum))/distr;}//pp_sigma and pp_pi diagonal
  else if((a!=b) && (a*b!=0)){h=(V[2]-V[3])*(vconum[a-1]*(*d).at(b-1)+vconum[b-1]*(*d).at(a-1)-2*(*d).at(a-1)*(*d).at(b-1)*(*d).at(conum))/distr;}//pp_sigma and pp_pi off-diagonal

  //V&G routine ends
  return h;
} //Hamder() ends

#endif

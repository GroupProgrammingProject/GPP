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
#include "functions.h"
#include "Gethijab.h"
#include <vector>

//Vector syntax: std::vector<double>. Assignments: xi=(*v).at(i). Inputs as pointer: std::vector<double>*

/*Change vector syntax as: std::vector<double>. Assignments xi=(*v).at(i). Inputs as pointer: std::vector<double>*, to call in functs ad &*/

double verlet(int norbs, double m, double rc,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* vx, std::vector<double>*vy, std::vector<double>* vz, Eigen::MatrixXi* c,double sx,double sy,double sz,std::vector<double>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXi* rx, Eigen::MatrixXi* ry, Eigen::MatrixXi* rz, Eigen::MatrixXi* modr);/*Inputs, in order: #orbitals; mass; cut-off radius; timestep; xyz arrays; velocity arrays; matrix of eigenvectors (N*N) (vectors as columns); cell sizes in xyz (put big numbers if you don't want PBCs); nearest neighbour lists; atom vector distances; modulus of vector distances.*/

void forces(int N,int norbs, Eigen::MatrixXi* rx, Eigen::MatrixXi* ry, Eigen::MatrixXi* rz, Eigen::MatrixXi* modr, Eigen::MatrixXi* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz);/*Inputs, in order: #atoms; #orbitals; atom vector distances; modulus of vector distances; matrix of eigenvectors (N*N) (vectors as columns); nearest neighbour lists; forces vectors.*/

void velocity(int N, double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T);/*Inputs, in order: #atoms; mass; vectors of velocities; temperature.*/

double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr,int conum);/*Inputs, in order: first atom index; second atom index; first orbital index; second orbital index; vector of direction cosines; modulus of distance between i and j; switch for component.*/ 

double verlet(int norbs, double m, double rc,double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, Eigen::MatrixXi* c,double sx,double sy,double sz,std::vector<double>* nnear,Eigen::MatrixXi* inear, Eigen::MatrixXi* rx, Eigen::MatrixXi* ry, Eigen::MatrixXi* rz, Eigen::MatrixXi* modr)
{ double boltz=1./11603;
  int N=(*x).size();
  std::vector<double> fx(N),fy(N),fz(N),fxn(N),fyn(N),fzn(N),vxm(N),vym(N),vzm(N);
  forces(N,norbs,rx,ry,rz,modr,c,nnear,inear,&fx,&fy,&fz); //calculate the forces
  
  for(int i=0; i<N; i++)
    {
      (*x).at(i)=(*x).at(i)+(*vx).at(i)*dt+0.5*fx.at(i)*dt*dt/m;
      (*y).at(i)=(*y).at(i)+(*vy).at(i)*dt+0.5*fy.at(i)*dt*dt/m;
      (*z).at(i)=(*z).at(i)+(*vz).at(i)*dt+0.5*fz.at(i)*dt*dt/m;
      forces(N,norbs,rx,ry,rz,modr,c,nnear,inear,&fx,&fy,&fz);//recalculate forces
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
      Tf=2*kin/(3*boltz*N); //final temperature
    }
  
  return Tf;	
}

/*INPUTS:N=numb. atoms; x,y,z=atom positions (arrays); c=eigenvectors (matrix with each column as the n-th eigenvector); rc=cut-off radius;
nnear=number of nearest neighbours (nn) to i-th atom (array); inear=label of j-th nn to i-th atom (matrix); fx,fy,fz forces on each atom (arrays);
maxnn=max number of nn */
void forces(int N,int norbs, Eigen::MatrixXi* rx, Eigen::MatrixXi* ry, Eigen::MatrixXi* rz, Eigen::MatrixXi* modr, Eigen::MatrixXi* c, std::vector<int>* nnear, Eigen::MatrixXi* inear, std::vector<double>* fx, std::vector<double>* fy, std::vector<double>* fz)
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
      sumphi=sumphi+o(r);
    }
    for(j=0;j<(*nnear).at(i);j++){ /*Cycle spanning the nearest neighbours of i*/
      nearlabel=(*inear)(i,j);
      for(m=0;m<(*nnear).at(nearlabel);m++){     
	dx=(*rx)(nearlabel,(*inear)(nearlabel,m)); /*Definition of vector distances*/
	dy=(*ry)(nearlabel,(*inear)(nearlabel,m));
	dz=(*rz)(nearlabel,(*inear)(nearlabel,m));
	r=(*modr)(nearlabel,(*inear)(nearlabel,m)); /*Modulus of distance*/
 
	sumphinn=sumphinn+o(r);
      }
      dx=(*rx)(i,nearlabel); /*Definition of vector distances*/
      dy=(*ry)(i,nearlabel);
      dz=(*rz)(i,nearlabel);
      
      r=(*modr)(i,nearlabel); /*Modulus of distance*/

      ddnorm.at(0)=dx/r;
      ddnorm.at(1)=dy/r;
      ddnorm.at(2)=dz/r;
      
      for(l=0;l<norbs;l++){ /*Cycle spanning the first orbital type*/
	for(lp=0;lp<norbs;lp++){ /*Cycle spanning the second orbital type*/
	  derivx=ds(r,dx)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,0);
	  derivy=ds(r,dy)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,1);
	  derivz=ds(r,dz)*Gethijab(i,nearlabel,l,lp,&ddnorm)+s(r)*Hamder(i,nearlabel,l,lp,&ddnorm,r,2);
	  for(n=0;n<norbs*N;n++){ /*Cycle spanning the level of the eigenvector*/
	    dualeigen=(*c).(n,l+i*norbs)*(*c)(n,lp+nearlabel*norbs);	
	    (*fx).at(i)=(*fx).at(i)-2*derivx*dualeigen;
	    (*fy).at(i)=(*fy).at(i)-2*derivy*dualeigen;
	    (*fz).at(i)=(*fz).at(i)-2*derivz*dualeigen;
	  }
	}
      }
      sumphinn=0;
      //calculation of repulsive forces
      (*fx).at(i)=(*fx).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dx);
      (*fy).at(i)=(*fy).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dy);
      (*fz).at(i)=(*fz).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*d_o(r,dz);
    }
    std::cout << (*fx).at(i) << std::endl;
    std::cout << (*fy).at(i) << std::endl;
    std::cout << (*fz).at(i) << std::endl;
  }  
}

void velocity(int N, double m, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T)
{  //calculate velocities
	double c,ke,boltz=1./11603;//1.38*pow(10,-23);
	double vxtot=0,vytot=0,vztot=0,msvx=0,msvy=0,msvz=0,vxavg,vyavg,vzavg,Tp;
	c=sqrt(3*boltz*T/m);
	for (int i=0; i<N; i++)
	{  
	  srand (1); //CHANGE RNG!!!!!!!!
	  (*vx).at(i)=c*(2*rand()-1);
	  (*vy).at(i)=c*(2*rand()-1);
	  (*vz).at(i)=c*(2*rand()-1);
	  vxtot=vxtot+(*vx).at(i);
	  vytot=vytot+(*vy).at(i);
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
	Tp=2*ke/(3*boltz*N);
	for(int i=0; i<N; i++)
	{
		(*vx).at(i)=(*vx).at(i)*sqrt(T/Tp);
		(*vy).at(i)=(*vy).at(i)*sqrt(T/Tp);
		(*vz).at(i)=(*vz).at(i)*sqrt(T/Tp);
	}
}
*/
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

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

int verlet(int N, int nmd, int norbs, std::vector<double>* mass, double rc, double rv, double T, double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* c,double sx,double sy,double sz,int nnmax);/*Inputs, in order: #atoms; #sim. steps; #orbitals; array containing masses (if all are equal
just define an array of dim 1); cut-off radius; verlet radius; initial temperature; timestep; xyz arrays; matrix of eigenvectors (N*N) (vectors as columns);
cell sizes in xyz (put big numbers if you don't want PBCs) */

void forces(int dim,int norbs,std::vector<double>* x,std::vector<double>* y,std::vector<double>* z,std::vector<double>* c,double rc,std::vector<int>* nnear,std::vector<int>* inear,std::vector<double>* fx,std::vector<double>* fy,std::vector<double>* fz);

void velocity(int N, std::vector<double>* mass, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T);

void near_neigh(int N, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, double rc, std::vector<int> *nnear, std::vector<int> *inear, double sx, double sy, double sz);



//using namespace std;
//
int verlet(int N, int nmd, int norbs, std::vector<double>* mass, double rc, double rv, double T, double dt, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,std::vector<double>* c,double sx,double sy,double sz,int nnmax)
{
	//Implement Velocity Verlet algorithm
//  double *fx, *fy, *fz, *vx, *vy, *vz, *xold, *yold, *zold, m=mass[1],*fxn,*fyn,*fzn,kin,*vxm,*vym,*vzm,dx,dy,dz,dist,dmax,svxm,svym,svzm,Tf;
		  std::vector<double> fx(N),fy(N),fz(N),vx(N),vy(N),vz(N),xold(N),yold(N),zold(N),fxn(N),fyn(N),fzn(N),vxm(N),vym(N),vzm(N);
		  double m=(*mass).at(0),kin,dx,dy,dz,dist,dmax,svxm,svym,svzm,Tf;
		 // int *nnear,*inear;
	std::vector<int> nnear(N), inear(N*nnmax);
  double boltz=1/11603;//Boltzmann's constant in eV/K  1.38*pow(10,-23);
//  fx= new double [N];
//	fy=new double [N];
//	fz=new double [N];
//	vx=new double [N];
//	vy=new double [N];
//	vz=new double [N];
//	xold=new double [N];
//	yold=new double [N];
//	zold=new double [N];
//	fxn=new double [N];
//	fyn=new double [N];
//	fzn=new double [N];
//	vxm=new double [N];
//	vym=new double [N];
//	vzm=new double [N];
//	nnear=new int [N];
//	inear=new int [N*nnmax];
	near_neigh(N,x,y,z,rc,&nnear,&inear,sx,sy,sz); 
	for(int i=0; i<N; i++)
	{
		xold.at(i)=(*x).at(i);
		yold.at(i)=(*y).at(i);//[i]=y[i];
		zold.at(i)=(*z).at(i);//zold[i]=z[i];
	}
	forces(N,norbs,x,y,z,c,rc,&nnear,&inear,&fx,&fy,&fz); //calculate the forces
	velocity(N,mass,&vx,&vy,&vz,T);
	
	for(int imd=1; imd<=nmd; imd++) //cycle through nmd steps
	{
		for(int i=0; i<N; i++)
		{
			(*x).at(i)=(*x).at(i)+vx.at(i)*dt+0.5*fx.at(i)*dt*dt/m;
			(*y).at(i)=(*y).at(i)+vy.at(i)*dt+0.5*fy.at(i)*dt*dt/m;//[i]=y[i]+vy[i]*dt+0.5*fy[i]*dt*dt/m;
			(*z).at(i)=(*z).at(i)+vz.at(i)*dt+0.5*fz.at(i)*dt*dt/m;//z[i]=z[i]+vz[i]*dt+0.5*fz[i]*dt*dt/m;

			dx=(*x).at(i)-xold.at(i);
			//dx=x[i]-xold[i];
			dx=dx-sx*round(dx/sx);
			dy=(*y).at(i)-yold.at(i);
			dy=dy-sy*round(dy/sy);
			dz=(*z).at(i)-zold.at(i);
			dz=dz-sz*round(dz/sz);
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			if (dist>dmax){dmax=dist;} //set dmax to largest dist so far
		}
		if(dmax>=0.4*(rv-rc)) //do the Verlet cages
		{
			near_neigh(N,x,y,z,rc,&nnear,&inear,sx,sy,sz);
			for(int i=0; i<N; i++)
			{
				xold.at(i)=(*x).at(i);
				yold.at(i)=(*y).at(i);
				zold.at(i)=(*z).at(i);
				//xold[i]=x[i];
				//yold[i]=y[i];
				//zold[i]=z[i];
			}
			forces(N,norbs,x,y,z,c,rc,&nnear,&inear,&fx,&fy,&fz);//recalculate forces
			for(int i=0; i<N; i++)//calculate new velocities
			{
				vx.at(i)=vx.at(i)+dt*(fx.at(i)+fxn.at(i))/(2*m);
				vy.at(i)=vy.at(i)+dt*(fy.at(i)+fyn.at(i))/(2*m);
				vz.at(i)=vz.at(i)+dt*(fz.at(i)+fzn.at(i))/(2*m);
//				vx[i]=vx[i]+dt*(fx[i]+fxn[i])/(2*m);
//				vy[i]=vy[i]+dt*(fy[i]+fyn[i])/(2*m);
//				vz[i]=vz[i]+dt*(fz[i]+fzn[i])/(2*m);
				fx.at(i)=fxn.at(i);
				fy.at(i)=fyn.at(i);
				fz.at(i)=fzn.at(i);
//				fx[i]=fxn[i];
//				fy[i]=fyn[i];
//				fz[i]=fzn[i];
			}
			for(int i=0; i<N; i++)//mean square velocities
			{
				svxm=svxm+vx.at(i)*vx.at(i);
				svym=svym+vy.at(i)*vy.at(i);
				svzm=svzm+vz.at(i)*vz.at(i);
				//svxm=svxm+vx[i]*vx[i];
				//svym=svym+vy[i]*vy[i];
				//svzm=svzm+vz[i]*vz[i];
			}
			kin=0,5*m*(svxm+svym+svzm); //kinetic energy
			Tf=2*kin/(3*boltz*N); //final temperature
		}
	}
}

/*INPUTS:dim=numb. atoms; x,y,z=atom positions (arrays); c=eigenvectors (matrix with each column as the n-th eigenvector); rc=cut-off radius;
nnear=number of nearest neighbours (nn) to i-th atom (array); inear=label of j-th nn to i-th atom (matrix); fx,fy,fz forces on each atom (arrays);
maxnn=max number of nn */
void forces(int dim,int norbs,std::vector<double>* x,std::vector<double>* y,std::vector<double>* z,std::vector<double>* c,double rc,std::vector<int>* nnear,std::vector<int>* inear,std::vector<double>* fx,std::vector<double>* fy,std::vector<double>* fz)
{ int k,i,j,l,lp,n,m,nearlabel; /* dummy indeces for cycles*/
  std::vector<double> dd(3),ddrx(3),ddrrx(3),ddlx(3),ddllx(3),ddry(3),ddrry(3),ddly(3),ddlly(3),ddrz(3),ddrrz(3),ddlz(3),ddllz(3);
//   double dd[3],ddrx[3],ddrrx[3],ddlx[3],ddllx[3],ddry[3],ddrry[3],ddrry[3],ddly[3],ddlly[3],ddrz[3],ddrrz[3],ddlz[3],ddllz[3]; 
  double ddm,ddmrx,ddmrrx,ddmlx,ddmllx,ddmry,ddmrry,ddmly,ddmlly,ddmrz,ddmrrz,ddmlz,ddmllz,h=0.001,sumphinn,sumphi,dualeigen,derivx,derivy,derivz;
	

  for(i=0;i<dim;i++){ /*initialisation of forces*/
    (*fx).at(i)=0;
    (*fy).at(i)=0;
    (*fz).at(i)=0;
    sumphinn=0;
  }

  

  for(i=0;i<dim;i++){ /*Cycle to compute band structure forces on atom i*/
    sumphi=0;
    for(k=0;k<(*nnear).at(i);k++){
      nearlabel=(*inear).at(i*10+k);
      dd.at(0)=(*x).at(i)-(*x).at(nearlabel); /*Definition of vector distances*/
      dd.at(1)=(*y).at(i)-(*y).at(nearlabel);
      dd.at(2)=(*z).at(i)-(*z).at(nearlabel);
      
      ddm=sqrt(dd.at(0)*dd.at(0)+dd.at(1)*dd.at(1)+dd.at(2)*dd.at(2)); /*Modulus of distance*/
      
      sumphi=sumphi+o(ddm);
   
    }
    for(j=0;j<(*nnear).at(i);j++){ /*Cycle spanning the nearest neighbours of i*/
      nearlabel=(*inear).at(i*10+j);
      dd.at(0)=(*x).at(i)-(*x).at(nearlabel); /*Definition of vector distances*/
	dd.at(1)=(*y).at(i)-(*y).at(nearlabel);
	dd.at(2)=(*z).at(i)-(*z).at(nearlabel);
	
	ddm=sqrt(dd.at(0)*dd.at(0)+dd.at(1)*dd.at(1)+dd.at(2)*dd.at(2)); /*Modulus of distance*/
	
	for(k=0;k<3;k++){ /*Initialisation of vector distances to perform derivatives */
	  ddrx.at(k)=dd.at(k); 
	  ddrrx.at(k)=dd.at(k);
	  ddlx.at(k)=dd.at(k);
	  ddllx.at(k)=dd.at(k);

	  ddry.at(k)=dd.at(k);
	  ddrry.at(k)=dd.at(k);
	  ddly.at(k)=dd.at(k);
	  ddlly.at(k)=dd.at(k);
	
	  ddrz.at(k)=dd.at(k);
	  ddrrz.at(k)=dd.at(k);
	  ddlz.at(k)=dd.at(k);
	  ddllz.at(k)=dd.at(k);
	}
	
	/*Definition of variables to take derivatives (r->x+h,rr->x+2h,l->x-h,ll->x-h)*/
	ddrx.at(0)=dd.at(0)+h;  
	ddrrx.at(0)=dd.at(0)+2*h;
	ddlx.at(0)=dd.at(0)-h;
	ddllx.at(0)=dd.at(0)-2*h;

	ddry.at(1)=dd.at(1)+h;  
	ddrry.at(1)=dd.at(1)+2*h;
	ddly.at(1)=dd.at(1)-h;
	ddlly.at(1)=dd.at(1)-2*h;

	ddrz.at(2)=dd.at(2)+h;  
	ddrrz.at(2)=dd.at(2)+2*h;
	ddlz.at(2)=dd.at(2)-h;
	ddllz.at(2)=dd.at(2)-2*h;

	ddmrx=sqrt(ddrx.at(0)*ddrx.at(0)+ddrx.at(1)*ddrx.at(1)+ddrx.at(2)*ddrx.at(2));
	ddmry=sqrt(ddry.at(0)*ddry.at(0)+ddry.at(1)*ddry.at(1)+ddry.at(2)*ddry.at(2));
	ddmrz=sqrt(ddrz.at(0)*ddrz.at(0)+ddrz.at(1)*ddrz.at(1)+ddrz.at(2)*ddrz.at(2));

	ddmrrx=sqrt(ddrrx.at(0)*ddrrx.at(0)+ddrrx.at(1)*ddrrx.at(1)+ddrrx.at(2)*ddrrx.at(2));
	ddmrry=sqrt(ddrry.at(0)*ddrry.at(0)+ddrry.at(1)*ddrry.at(1)+ddrry.at(2)*ddrry.at(2));
	ddmrrz=sqrt(ddrrz.at(0)*ddrrz.at(0)+ddrrz.at(1)*ddrrz.at(1)+ddrrz.at(2)*ddrrz.at(2));
      
	ddmlx=sqrt(ddlx.at(0)*ddlx.at(0)+ddlx.at(1)*ddlx.at(1)+ddlx.at(2)*ddlx.at(2));
	ddmly=sqrt(ddly.at(0)*ddly.at(0)+ddly.at(1)*ddly.at(1)+ddly.at(2)*ddly.at(2));
	ddmlz=sqrt(ddlz.at(0)*ddlz.at(0)+ddlz.at(1)*ddlz.at(1)+ddlz.at(2)*ddlz.at(2));

	ddmllx=sqrt(ddllx.at(0)*ddllx.at(0)+ddllx.at(1)*ddllx.at(1)+ddllx.at(2)*ddllx.at(2));
	ddmlly=sqrt(ddlly.at(0)*ddlly.at(0)+ddlly.at(1)*ddlly.at(1)+ddlly.at(2)*ddlly.at(2));
	ddmllz=sqrt(ddllz.at(0)*ddllz.at(0)+ddllz.at(1)*ddllz.at(1)+ddllz.at(2)*ddllz.at(2));

	for(k=0;k<3;k++){
	  ddrx.at(k)=ddrx.at(k)/ddmrx;
	  ddrrx.at(k)=ddrrx.at(k)/ddmrrx;
	  ddlx.at(k)=ddlx.at(k)/ddmlx;
	  ddllx.at(k)=ddllx.at(k)/ddmllx;

	  ddry.at(k)=ddry.at(k)/ddmry;
	  ddrry.at(k)=ddrry.at(k)/ddmrry;
	  ddly.at(k)=ddly.at(k)/ddmly;
	  ddlly.at(k)=ddlly.at(k)/ddmlly;

	  ddrz.at(k)=ddrz.at(k)/ddmrz;
	  ddrrz.at(k)=ddrrz.at(k)/ddmrrz;
	  ddlz.at(k)=ddlz.at(k)/ddmlz;
	  ddllz.at(k)=ddllz.at(k)/ddmllz;
	}

//	double h=Gethijab(i,j,l,lp,&ddrx,6,6);

	for(l=0;l<norbs;l++){ /*Cycle spanning the first orbital type*/
	  for(lp=0;lp<norbs;lp++){ /*Cycle spanning the second orbital type*/
	    //derivx=(-Gethijab(i,nearlabel,l,lp,&ddrrx)*s(ddmrrx)+8*Gethijab(i,nearlabel,l,lp,&ddrx)*s(ddmrx)-8*Gethijab(i,nearlabel,l,lp,&ddlx)*s(ddmlx)+Gethijab(i,nearlabel,l,lp,&ddllx)*s(ddmllx))/(12*h);
	      //derivy=(-Gethijab(i,nearlabel,l,lp,&ddrry)*s(ddmrry)+8*Gethijab(i,nearlabel,l,lp,&ddry)*s(ddmry)-8*Gethijab(i,nearlabel,l,lp,&ddly)*s(ddmly)+Gethijab(i,nearlabel,l,lp,&ddlly)*s(ddmlly))/(12*h);
	      //derivz=(-Gethijab(i,nearlabel,l,lp,&ddrrz)*s(ddmrrz)+8*Gethijab(i,nearlabel,l,lp,&ddrz)*s(ddmrz)-8*Gethijab(i,nearlabel,l,lp,&ddlz)*s(ddmlz)+Gethijab(i,nearlabel,l,lp,&ddllz)*s(ddmllz))/(12*h);

	    derivx=(Gethijab(i,nearlabel,l,lp,&ddrx)*s(ddmrx)-Gethijab(i,nearlabel,l,lp,&ddlx)*s(ddmlx))/(2*h);
	    derivy=(Gethijab(i,nearlabel,l,lp,&ddry)*s(ddmry)-Gethijab(i,nearlabel,l,lp,&ddly)*s(ddmly))/(2*h);
	    derivz=(Gethijab(i,nearlabel,l,lp,&ddrz)*s(ddmrz)-Gethijab(i,nearlabel,l,lp,&ddlz)*s(ddmlz))/(2*h);

	    for(n=0;n<dim;n++){ /*Cycle spanning the level of the eigenvector*/
		 dualeigen=(*c).at(l+i+n*4*dim)*(*c).at(lp+nearlabel+n*4*dim);

		 (*fx).at(i)=(*fx).at(i)-2*2*derivx*dualeigen;
		 (*fy).at(i)=(*fy).at(i)-2*2*derivy*dualeigen;
		 (*fz).at(i)=(*fz).at(i)-2*2*derivz*dualeigen;
				
	       }
	  }
	}


	std::cout << "fx(" << i << ")=" << std::setprecision(10) << (*fx).at(i) << std::endl;
	std::cout << "fy(" << i << ")=" << std::setprecision(10) << (*fy).at(i) << std::endl;
	std::cout << "fz(" << i << ")=" << std::setprecision(10) << (*fz).at(i) << std::endl;
	std::cout << std::endl;

	sumphinn=0;

	for(m=0;m<(*nnear).at(nearlabel);m++){
	  dd.at(0)=(*x).at(nearlabel)-(*x).at((*inear).at(nearlabel*dim+m)); /*Definition of vector distances*/
	  dd.at(1)=(*y).at(nearlabel)-(*y).at((*inear).at(nearlabel*dim+m));
	  dd.at(2)=(*z).at(nearlabel)-(*z).at((*inear).at(nearlabel*dim+m));
		    
	    ddm=sqrt(dd.at(0)*dd.at(0)+dd.at(1)*dd.at(1)+dd.at(2)*dd.at(2)); /*Modulus of distance*/
	    
	    sumphinn=sumphinn+o(ddm);
	}
	      /*calculation of repuslve forces*/
	(*fx).at(i)=(*fx).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrrx)+8*o(ddmrx)-8*o(ddmlx)+o(ddmllx))/(12*h);
	(*fy).at(i)=(*fy).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrry)+8*o(ddmry)-8*o(ddmly)+o(ddmlly))/(12*h);
	(*fz).at(i)=(*fz).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrrz)+8*o(ddmrz)-8*o(ddmlz)+o(ddmllz))/(12*h);
	
      
    }
  }
}

void near_neigh(int N, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, double rc, std::vector<int> *nnear, std::vector<int> *inear, double sx, double sy, double sz)
{  //determine the nearest neighbours for each atom
	double dx,dy,dz,dist;
	for (int i=0; i<N; i++) {(*nnear).at(i)=0;}// nnear[i]=0; }
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			dx=(*x).at(i)-(*x).at(j);
			dy=(*y).at(i)-(*y).at(j);
			dz=(*z).at(i)-(*z).at(j);
			//dx=x[i]-x[j];
			dx=dx-sx*round(dx/sx);//if dx>sx, the distance to N.N. is dx-sx
//			std::cout << "check that new dx isn't zero, dx= " << dx << std::endl;
			//dy=y[i]-y[j];
			dy=dy-sy*round(dy/sy);
			//dz=z[i]-z[j];
			dz=dz-sz*round(dz/sz);
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			if (dist<rc && i!=j) //add the atom j to the nearest neighbour list of i if this holds
			{
			(*nnear).at(i)++;
			(*inear).at(i*10+(*nnear).at(i)-1)=j;;
			//	nnear[i]=nnear[i]+1;
			//	inear[i*N+nnear[i]]=j; //a matrix with i rows, nnear[i] (no of nearest neighbours) columns
			}
		}
	}
}	//ends near_neigh

void velocity(int N, std::vector<double>* mass, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T)
{  //calculate velocities
	double c, ke, m=(*mass).at(1),boltz=1/11603;//1.38*pow(10,-23); //for now, each mass is the same
	double vxtot=0,vytot=0,vztot=0,msvx=0,msvy=0,msvz=0,vxavg,vyavg,vzavg,Tp;
	double g=sqrt(3*boltz*T/m);
	for (int i=0; i<N; i++)
	{
	  
	  srand (1);
	  (*vx).at(i)=0; (*vy).at(i)=0; (*vz).at(i)=0;
	  (*vx).at(i)=c*(2*rand()-1);
	  (*vy).at(i)=c*(2*rand()-1);
	  (*vz).at(i)=c*(2*rand()-1);
	  vxtot=vxtot+(*vx).at(i);
	  vytot=vytot+(*vy).at(i);
//		vx[i]=0; vy[i]=0; vz[i]=0;
//		c[i]=sqrt(3*boltz*T/m[i])
//		vx[i]=c[i]*(2*rand(0)-1); //random number generator required
//		vy[i]=c[i]*(2*rand(0)-1);
//		vz[i]=c[i]*(2*rand(0)-1);
//		vx[i]=c*(2*rand()-1);
//		vy[i]=c*(2*rand()-1);
//		vz[i]=c*(2*rand()-1);
//		vxtot=vxtot+vx[i];
//		vytot=vytot+vy[i];
//		vztot=vztot+vz[i];
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
//		vx[i]=vx[i]-vxavg;
//		vy[i]=vy[i]-vyavg;
//		vz[i]=vz[i]-vzavg;
//		msvx=msvx+vx[i]*vx[i];//mean square velocity
//		msvy=msvy+vy[i]*vy[i];
//		msvz=msvz+vz[i]*vz[i];
	}
	ke=0.5*m*(msvx+msvy+msvz); //must be adapted for different masses
	Tp=2*ke/(3*boltz);
	for(int i=0; i<N; i++)
	{
		(*vx).at(i)=(*vx).at(i)*sqrt(T/Tp);
		(*vy).at(i)=(*vy).at(i)*sqrt(T/Tp);
		(*vz).at(i)=(*vz).at(i)*sqrt(T/Tp);

//		vx[i]=vx[i]*sqrt(T/Tp);
//		vy[i]=vy[i]*sqrt(T/Tp);
//		vz[i]=vz[i]*sqrt(T/Tp);
	}
}

/*
       program dinamica molecolare
	implicit real*8 (a-h,o-z) 
        parameter (ntot=432)
	parameter (nvicmax=200)
	parameter (rc=4.5d0) ! cut-off radius 
	parameter (rv=5.d0)
	parameter (epsi=0.345)
	parameter (sigma=2.644)
	parameter (sx=17.64739)
	parameter (sy=17.64739)
	parameter (sz=24.957179*10)
	parameter (T=584.d0) 
	parameter (nmd=4000)
	parameter (dt=1.d-15) !timestep
	dimension x(ntot),y(ntot),z(ntot),nvicini(ntot) !nvicini(i)-> number of n.n. to i-th atom
	dimension ivicini(ntot,nvicmax) !ivicini(i,j)-> j-th n.n. to i-th atom 
	dimension parz(ntot)
	dimension fx(ntot),fy(ntot),fz(ntot)
	dimension fxn(ntot),fyn(ntot),fzn(ntot)
	dimension vx(ntot),vy(ntot),vz(ntot)
	dimension xsave(ntot),ysave(ntot),zsave(ntot)
	real*4 elapsed(2)
	real*4 total
        
	iseed=9
	aa=rand(iseed)
   	
	xk=1./11603.            !Boltzmann constant
	xm=108.*(1.66d-27)/16.	!Atomic mass

	rprimo=rc-0.3d0		!Polynomial raccord to potential
        diffe=rprimo-rc
	plus=rprimo+rc
        Ee=4.*epsi*((sigma/rprimo)**12.-(sigma/rprimo)**6.) 
	F=24.*epsi*((sigma**6.)/(rprimo**7.)-(2.*(sigma**12.)/
     +  (rprimo**13.)))
	D=F/(diffe**2.)-(2.*Ee/(diffe**3.))
	C=F/(2.*diffe)-(1.5d0)*D*plus
	B=-2.*C*rc-3.*D*(rc**2.)
	A=C*(rc**2.)+2.*D*(rc**3.)

	iprogstep=nmd/10
	ipercent=10     

	call leggipos(ntot,x,y,z)  !call to the subroutine reading atomic positions
	call vicini(ntot,nvicmax,x,y,z,rv,nvicini,ivicini,sx,sy,sz) !call to the subroutine calculating distances and n.n.s	
	
	do i=1,ntot
	 xsave(i)=x(i)
	 ysave(i)=y(i)
	 zsave(i)=z(i)
	end do

 	call ene(x,y,z,nvicini,ivicini,Epot,parz,epsi,sigma,ntot,
     +  nvicmax,rprimo,rc,A,B,C,D,sx,sy,sz) !call to energy calculation subroutine 
	call forza(x,y,z,nvicini,ivicini,fx,fy,fz,epsi,sigma,ntot
     +   ,nvicmax,rprimo,rc,A,B,C,D,sx,sy,sz) !call to force calculation subroutine
        call velini (xk,T,xm,vx,vy,vz,ntot,Tprimo,vxm,vym,vzm)

	dmax=0.d0

        ncall=0

	do imd=1,nmd   
	 do i=1,ntot
	  x(i)=x(i)+vx(i)*dt+0.5*(fx(i)/xm)*dt**2 !velocity verlet algorithm
	  y(i)=y(i)+vy(i)*dt+0.5*(fy(i)/xm)*dt**2
	  z(i)=z(i)+vz(i)*dt+0.5*(fz(i)/xm)*dt**2

	  ddx=x(i)-xsave(i)
	  ddx=ddx-sx*anint(ddx/sx)
	  ddy=y(i)-ysave(i)
	  ddy=ddy-sy*anint(ddy/sy)
	  ddz=z(i)-zsave(i)
	  ddz=ddz-sz*anint(ddz/sz)
          dd=sqrt(ddx*ddx+ddy*ddy+ddz*ddz)

          if(dd.ge.dmax)dmax=dd
	 end do
 
         if(dmax.ge.0.4*(rv-rc))then
          
          call vicini(ntot,nvicmax,x,y,z,rv,nvicini,ivicini,sx,sy,sz)
          ncall=ncall+1
          dmax=0.d0
         
          do l=1,ntot
           xsave(l)=x(l)
	   ysave(l)=y(l)
	   zsave(l)=z(l)
          end do

 	 end if
	 
	 call forza(x,y,z,nvicini,ivicini,fxn,fyn,fzn,epsi,sigma,ntot    !recalculation of forces
     +   ,nvicmax,rprimo,rc,A,B,C,D,sx,sy,sz)

         do i=1,ntot !cycle to recalculate atomic velocities
	  vx(i)=vx(i)+(0.5/xm)*(fxn(i)+fx(i))*dt
	  vy(i)=vy(i)+(0.5/xm)*(fyn(i)+fy(i))*dt
	  vz(i)=vz(i)+(0.5/xm)*(fzn(i)+fz(i))*dt
	  fx(i)=fxn(i)
	  fy(i)=fyn(i)
	  fz(i)=fzn(i)
	 end do
        
	call ene(x,y,z,nvicini,ivicini,Epot,parz,epsi,sigma,ntot,     !recalculation of potential energy
     +  nvicmax,rprimo,rc,A,B,C,D,sx,sy,sz)

        svxm=0
	svym=0
	svzm=0

        do i=1,ntot 
	 svxm=svxm+vx(i)**2
	 svym=svym+vy(i)**2
	 svzm=svzm+vz(i)**2
	end do

	EnCin=(0.5)*xm*(svxm+svym+svzm)

	Tfin=(2.*EnCin/ntot)/(3.*xk)

	 write (20,*) dt*imd,Tfin,EnCin,Epot,EnCin+Epot

	if(imd.ge.iprogstep)Then	
	iprogress=mod(imd,iprogstep)
         if(iprogress.eq.0)Then   
           print *,ipercent,'% complete'
	   ipercent=ipercent+10
         end if
	end if		

        end do

	

	print *,'Initial temperature:'
	print *, T,'K'
	print *, 'Timestep:'
	print *, dt,'secondi'
	print *, 'Number of steps'
	write (*,*) nmd
	print *, 'Total simulation time:'
	print *, nmd*dt,'secondi'
	
	print *,'t vs T, E cin, E pot, E tot -> .20'

	total=etime(elapsed)
	minutes=aint(total/60)
	seconds=total-minutes*60
	
	stop
	end

!*******************************************************!
!***  READING ATOMIC POSITIONS  ************************!
!*******************************************************!
	subroutine leggipos(ntot,x,y,z) 
	implicit real*8 (a-h,o-z)
	dimension x(ntot),y(ntot),z(ntot)
	open(30,file='fccseconditerzi.txt',status='old') 
	
	do i=1,ntot
	 read(30,*)x(i),y(i),z(i) 
	end do

	close(30)
        return
	end

!*******************************************************!
!***  CALCULATION OF N.N.S  ****************************!
!*******************************************************!
	subroutine vicini(ntot,nvicmax,x,y,z,rc,nvicini,ivicini,
     +  sx,sy,sz) 
	implicit real*8 (a-h,o-z)
	dimension x(ntot),y(ntot),z(ntot),nvicini(ntot)
	dimension ivicini(ntot,nvicmax)
        
	do i=1,ntot
         nvicini(i)=0
	end do
	
	do i=1,ntot
	 
	 do j=1,ntot 
	  dx=x(i)-x(j)
	  dx=dx-sx*anint(dx/sx) !P.B.C. implementation
	  dx2=dx*dx
	  dy=y(i)-y(j)
	  dy=dy-sy*anint(dy/sy)
	  dy2=dy*dy
	  dz=z(i)-z(j)
	  dz=dz-sz*anint(dz/sz)
	  dz2=dz*dz
	  dist2=dx2+dy2+dz2
	  dist=dsqrt(dist2) 
	  
	  if(dist.le.rc.and.i.ne.j)Then !.le.-> less or equal; ;.ne.-> not equal
           nvicini(i)=nvicini(i)+1 !Increment of the number of n.n.s
	   ivicini(i,nvicini(i))=j !labelling of n.n.s to the i-th atom
	  endif

       	 end do 
	
 	end do 		
	
	return
	end

!*******************************************************!
!***  CALCOLO ENERGIE  *********************************!
!*******************************************************!
	subroutine ene(x,y,z,nvicini,ivicini,Epot,parz,epsi,sigma,ntot !dichiaro la subroutine di calcolo dell'energia
     +	,nvicmax,rprimo,rc,A,B,C,D,sx,sy,sz)
	implicit real*8 (a-h,o-z)
	dimension x(ntot),y(ntot),z(ntot),nvicini(ntot)
	dimension ivicini(ntot,nvicmax)
	dimension parz(ntot)
        Epot=0	
	do i=1,ntot
	 parz(i)=0	
	 
	 do j=1,nvicini(i)
	  k=ivicini(i,j)
	  dx=x(i)-x(k)
	  dx=dx-sx*anint(dx/sx)
	  dx2=dx*dx
	  dy=y(i)-y(k)
	  dy=dy-sy*anint(dy/sy)
	  dy2=dy*dy
	  dz=z(i)-z(k)
	  dz=dz-sz*anint(dz/sz)
	  dz2=dz*dz
	  dist2=dx2+dy2+dz2
	  dist=dsqrt(dist2)
          
	  if(dist.le.rprimo)Then
	   parz(i)=parz(i)+4.*epsi*((sigma/dist)**12.-(sigma/dist)**6.)
	  else
            if(dist.le.rc.and.dist.gt.rprimo)Then
	     parz(i)=parz(i)+A+B*dist+C*(dist**2)+D*(dist**3) !Raccord polynomial
            else
	     continue
	    end if
	  end if
	  
	  end do

	 Epot=Epot+parz(i)
	
	end do
	
	Epot=Epot/2

	return
	end

!*******************************************************!
!***  CALCULATION OF FORCES  ***************************!
!*******************************************************!
	subroutine forza(x,y,z,nvicini,ivicini,fx,fy,fz,epsi,sigma,ntot !subroutine di calcolo forze
     +  ,nvicmax,rprimo,rc,A,B,C,D,sx,sy,sz)
        implicit real*8 (a-h,o-z)
	dimension x(ntot),y(ntot),z(ntot),nvicini(ntot)
	dimension ivicini(ntot,nvicmax)
        dimension fx(ntot),fy(ntot),fz(ntot)

	do i=1,ntot
	 fx(i)=0
	 fy(i)=0
	 fz(i)=0

	 do j=1,nvicini(i)
          k=ivicini(i,j)
	  dx=x(i)-x(k)
	  dx=dx-sx*anint(dx/sx)
	  dx2=dx*dx
	  dy=y(i)-y(k)
	  dy=dy-sy*anint(dy/sy)
	  dy2=dy*dy
	  dz=z(i)-z(k)
	  dz=dz-sz*anint(dz/sz)
	  dz2=dz*dz
	  dist2=dx2+dy2+dz2
	  dist=dsqrt(dist2)

	  if(dist.le.rprimo)Then
	  fx(i)=fx(i)+24*epsi*((sigma**6)/(dist**8))*(2*((sigma/dist)**6
     +    )-1)*dx
          fy(i)=fy(i)+24*epsi*((sigma**6)/(dist**8))*(2*((sigma/dist)**6
     +    )-1)*dy
	  fz(i)=fz(i)+24*epsi*((sigma**6)/(dist**8))*(2*((sigma/dist)**6
     +    )-1)*dz
	  
	  else
            if(dist.le.rc.and.dist.gt.rprimo)Then
	     fx(i)=fx(i)-dx*((B+3*D*dist2)/(dist)+2*C)
	     fy(i)=fy(i)-dy*((B+3*D*dist2)/(dist)+2*C)
	     fz(i)=fz(i)-dz*((B+3*D*dist2)/(dist)+2*C)
          
            else
	     continue
	    end if
	  end if	  
          
	 end do

	end do

        return
	end

!*******************************************************!
!***  CALCULATION OF INITIAL VELOCITIES  ***************!
!*******************************************************!
	subroutine velini (xk,T,xm,vx,vy,vz,ntot,Tprimo,vxm,vym,vzm)
	implicit real*8 (a-h,o-z)
	dimension vx(ntot),vy(ntot),vz(ntot)
 
	c=dsqrt(((3*xk*T)/xm))

	vxm=0
	vym=0
	vzm=0
        svxm=0
	svym=0
	svzm=0
	do i=1,ntot
	 vx(i)=0
	 vx(i)=c*(2*rand(0)-1)
	 vy(i)=0
 	 vy(i)=c*(2*rand(0)-1)
	 vz(i)=0
	 vz(i)=c*(2*rand(0)-1)
	 vxm=vxm+vx(i)
	 vym=vym+vy(i)
	 vzm=vzm+vz(i)
	end do

 	vxm=vxm/ntot
	vym=vym/ntot
	vzm=vzm/ntot

        do i=1,ntot	   !velocity rescaling to impose zero overall momentum
	 vx(i)=vx(i)-vxm
	 vy(i)=vy(i)-vym
	 vz(i)=vz(i)-vzm
	 svxm=svxm+vx(i)**2
	 svym=svym+vy(i)**2
	 svzm=svzm+vz(i)**2
	end do

	svxm=svxm/ntot
	svym=svym/ntot
	svzm=svzm/ntot

        EnCin=0.5*xm*(svxm+svym+svzm)

	Tprimo=(2*EnCin)/(3*xk)

	svxm=0
	svym=0
	svzm=0

	do i=1,ntot	   
	 vx(i)=vx(i)*dsqrt(T/Tprimo)
	 vy(i)=vy(i)*dsqrt(T/Tprimo)
	 vz(i)=vz(i)*dsqrt(T/Tprimo)
	end do  !rescaling imposing a definite temperature

	return
  	end

*/

#endif

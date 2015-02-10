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
  double boltz=1./11603;//Boltzmann's constant in eV/K  1.38*pow(10,-23);
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
{ int k,i,j,l,lp,n,m; /* dummy indeces for cycles*/
  std::vector<double> dd(3),ddrx(3),ddrrx(3),ddlx(3),ddllx(3),ddry(3),ddrry(3),ddly(3),ddlly(3),ddrz(3),ddrrz(3),ddlz(3),ddllz(3);
//   double dd[3],ddrx[3],ddrrx[3],ddlx[3],ddllx[3],ddry[3],ddrry[3],ddrry[3],ddly[3],ddlly[3],ddrz[3],ddrrz[3],ddlz[3],ddllz[3]; 
	double ddm,ddmrx,ddmrrx,ddmlx,ddmllx,ddmry,ddmrry,ddmly,ddmlly,ddmrz,ddmrrz,ddmlz,ddmllz,h=rc/10000,sumphinn,sumphi;

  for(i=0;i<dim;i++){ /*initialisation of forces*/
    (*fx).at(i)=0;
    (*fy).at(i)=0;
    (*fz).at(i)=0;
  }

  for(i=0;i<dim;i++){ /*Cycle to compute band structure forces on atom i*/
    sumphi=0;
    for(k=0;k<(*nnear).at(i);k++){
      dd.at(0)=((*x).at(i)-(*x).at((*inear).at(i*dim+k))); /*Definition of vector distances*/
      dd.at(1)=((*y).at(i)-(*y).at((*inear).at(i*dim+k)));
      dd.at(2)=((*z).at(i)-(*z).at((*inear).at(i*dim+k)));
      
      
      ddm=sqrt(dd.at(0)*dd.at(0)+dd.at(1)*dd.at(1)+dd.at(2)*dd.at(2)); /*Modulus of distance*/
      
      sumphi=sumphi+o(ddm);
    }
    for(j=0;j<(*nnear).at(i);j++){ /*Cycle spanning the nearest neighbours of i*/
      /*Check to avoid self interaction (redundant, but saves operations)*/
	dd.at(0)=((*x).at(i)-(*x).at((*inear).at(i*dim+j))); /*Definition of vector distances*/
	dd.at(1)=((*y).at(i)-(*y).at((*inear).at(i*dim+j)));
	dd.at(2)=((*z).at(i)-(*z).at((*inear).at(i*dim+j)));

	
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
	ddrx.at(0)=ddrx.at(0)+h; 
	ddrrx.at(0)=ddrrx.at(0)+2*h;
	ddlx.at(0)=ddlx.at(0)-h;
	ddllx.at(0)=ddllx.at(0)-2*h;

	ddry.at(1)=ddry.at(1)+h; 
	ddrry.at(1)=ddrry.at(1)+2*h;
	ddly.at(1)=ddly.at(1)-h;
	ddlly.at(1)=ddlly.at(1)-2*h;

	ddrz.at(2)=ddrz.at(2)+h; 
	ddrrz.at(2)=ddrrz.at(2)+2*h;
	ddlz.at(2)=ddlz.at(2)-h;
	ddllz.at(2)=ddllz.at(2)-2*h;

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

	std::cout << "fine" << std::endl;

	for(l=0;l<norbs;l++){ /*Cycle spanning the first orbital type*/
	  for(lp=0;lp<norbs;l++){ /*Cycle spanning the second orbital type*/
	    for(n=0;n<(dim);n++){ /*Cycle spanning the level of the eigenvector*/
	      (*fx).at(i)=(*fx).at(i)-2*(-Gethijab(i,j,l,lp,&ddrrx,6,6)*s(ddmrrx,6,6)+8*Gethijab(i,j,l,lp,&ddrx,6,6)*s(ddmrx,6,6)-8*Gethijab(i,j,l,lp,&ddlx,6,6)*s(ddmlx,6,6)+Gethijab(i,j,l,lp,&ddllx,6,6)*s(ddmllx,6,6))/(12*h)*(*c).at((l+i)+n*dim)*(*c).at((lp+j)+n*dim);

	      (*fy).at(i)=(*fy).at(i)-2*(-Gethijab(i,j,l,lp,&ddrry,6,6)*s(ddmrry,6,6)+8*Gethijab(i,j,l,lp,&ddry,6,6)*s(ddmry,6,6)-8*Gethijab(i,j,l,lp,&ddly,6,6)*s(ddmly,6,6)+Gethijab(i,j,l,lp,&ddlly,6,6)*s(ddmlly,6,6))/(12*h)*(*c).at((l+i)+n*dim)*(*c).at((lp+j)+n*dim);

	      (*fz).at(i)=(*fz).at(i)-2*(-Gethijab(i,j,l,lp,&ddrrz,6,6)*s(ddmrrz,6,6)+8*Gethijab(i,j,l,lp,&ddrz,6,6)*s(ddmrz,6,6)-8*Gethijab(i,j,l,lp,&ddlz,6,6)*s(ddmlz,6,6)+Gethijab(i,j,l,lp,&ddllz,6,6)*s(ddmllz,6,6))/(12*h)*(*c).at((l+i)+n*dim)*(*c).at((lp+j)+n*dim);

	      	std::cout << i << " " << j << " " << l << " " << lp << ": " << Gethijab(i,j,l,lp,&dd,6,6) << std::endl;
		std::cout << j << " " << i << " " << l << " " << lp << ": " << Gethijab(j,i,l,lp,&dd,6,6) << std::endl;
	    }
	  }
	}


	

	//	std::cout << (*fx).at(i) << std::endl;
	//	std::cout << (*fy).at(i) << std::endl;
	//	std::cout << (*fz).at(i) << std::endl;
	
	
	sumphinn=0;

	for(m=0;m<(*nnear).at((*inear).at(i*dim+j));m++){
	  dd.at(0)=abs((*x).at((*inear).at(i*dim+j))-(*x).at((*inear).at((*inear).at(i*dim+j)*dim+m))); /*Definition of vector distances*/
	  dd.at(1)=abs((*y).at((*inear).at(i*dim+j))-(*y).at((*inear).at((*inear).at(i*dim+j)*dim+m)));
	  dd.at(2)=abs((*z).at((*inear).at(i*dim+j))-(*z).at((*inear).at((*inear).at(i*dim+j)*dim+m)));
		    
	    ddm=sqrt(dd.at(0)*dd.at(0)+dd.at(1)*dd.at(1)+dd.at(2)*dd.at(2)); /*Modulus of distance*/
	    
	    sumphinn=sumphinn+o(ddm);
	}
	      /*calculation of repuslve forces*/
	(*fx).at(i)=(*fx).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrrx)+8*o(ddmrx)-8*o(ddmlx)+o(ddmllx))/(12*h);
	(*fy).at(i)=(*fy).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrry)+8*o(ddmry)-8*o(ddmly)+o(ddmlly))/(12*h);
	(*fz).at(i)=(*fz).at(i)-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrrz)+8*o(ddmrz)-8*o(ddmlz)+o(ddmllz))/(12*h);
	
	
      
    }
  }

//  return 0;

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
			(*inear).at(i*N+(*nnear).at(i)-1)=j;
		     
			}
		}
	}

}	//ends near_neigh

void velocity(int N, std::vector<double>* mass, std::vector<double>* vx, std::vector<double>* vy, std::vector<double>* vz, double T)
{  //calculate velocities
	double c, ke, m=(*mass).at(1),boltz=1./11603;//1.38*pow(10,-23); //for now, each mass is the same
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



#endif

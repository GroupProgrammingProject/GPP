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
#include "functions.h"
#include "Gethijab.h"

int verlet(int N, int nmd, int norbs, double *mass, double rc, double rv, double T, double dt, double *x, double *y, double *z, double *c,double sx,double sy,double sz,int nnmax);/*Inputs, in order: #atoms; #sim. steps; #orbitals; array containing masses (if all are equal
just define an array of dim 1); cut-off radius; verlet radius; initial temperature; timestep; xyz arrays; matrix of eigenvectors (N*N) (vectors as columns);
cell sizes in xyz (put big numbers if you don't want PBCs) */

int forces(int dim,int norbs,double *x,double *y,double *z,double *c,double rc,int *nnear,int *inear,double *fx,double *fy,double *fz);

void velocity(int N, double *mass, double *vx, double *vy, double *vz, double T);

void near_neigh(int N, double *x, double *y, double *z, double rc, int *nnear, int *inear, double sx, double sy, double sz);



using namespace std;
//
int verlet(int N, int nmd, int norbs, double *mass, double rc, double rv, double T, double dt, double *x, double *y, double *z,double *c,double sx,double sy,double sz,int nnmax)
{
	//Implement Velocity Verlet algorithm
  double *fx, *fy, *fz, *vx, *vy, *vz, *xold, *yold, *zold, m=mass[1],*fxn,*fyn,*fzn,kin,*vxm,*vym,*vzm,dx,dy,dz,dist,dmax,svxm,svym,svzm,Tf;
  int *nnear,*inear;
  double boltz=1/11603;//Boltzmann's constant in eV/K  1.38*pow(10,-23);
  fx= new double [N];
	fy=new double [N];
	fz=new double [N];
	vx=new double [N];
	vy=new double [N];
	vz=new double [N];
	xold=new double [N];
	yold=new double [N];
	zold=new double [N];
	fxn=new double [N];
	fyn=new double [N];
	fzn=new double [N];
	vxm=new double [N];
	vym=new double [N];
	vzm=new double [N];
	nnear=new int [N];
	inear=new int [N*nnmax];
	near_neigh(N,x,y,z,rc,nnear,inear,sx,sy,sz); 
	for(int i=0; i<N; i++)
	{
		xold[i]=x[i];
		yold[i]=y[i];
		zold[i]=z[i];
	}
	forces(N,norbs,x,y,z,c,rc,nnear,inear,fx,fy,fz); //calculate the forces
	velocity(N,mass,vx,vy,vz,T);
	
	for(int imd=1; imd<=nmd; imd++) //cycle through nmd steps
	{
		for(int i=0; i<N; i++)
		{
			x[i]=x[i]+vx[i]*dt+0.5*fx[i]*dt*dt/m;
			y[i]=y[i]+vy[i]*dt+0.5*fy[i]*dt*dt/m;
			z[i]=z[i]+vz[i]*dt+0.5*fz[i]*dt*dt/m;

			dx=x[i]-xold[i];
			dx=dx-sx*round(dx/sx);
			dy=y[i]-yold[i];
			dy=dy-sy*round(dy/sy);
			dz=z[i]-zold[i];
			dz=dz-sz*round(dz/sz);
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			if (dist>dmax){dmax=dist;} //set dmax to largest dist so far
		}
		if(dmax>=0.4*(rv-rc)) //do the Verlet cages
		{
			near_neigh(N,x,y,z,rc,nnear,inear,sx,sy,sz);
			for(int i=0; i<N; i++)
			{
				xold[i]=x[i];
				yold[i]=y[i];
				zold[i]=z[i];
			}
			forces(N,norbs,x,y,z,c,rc,nnear,inear,fx,fy,fz);//recalculate forces
			for(int i=0; i<N; i++)//calculate new velocities
			{
				vx[i]=vx[i]+dt*(fx[i]+fxn[i])/(2*m);
				vy[i]=vy[i]+dt*(fy[i]+fyn[i])/(2*m);
				vz[i]=vz[i]+dt*(fz[i]+fzn[i])/(2*m);
				fx[i]=fxn[i];
				fy[i]=fyn[i];
				fz[i]=fzn[i];
			}
			for(int i=0; i<N; i++)//mean square velocities
			{
				svxm=svxm+vx[i]*vx[i];
				svym=svym+vy[i]*vy[i];
				svzm=svzm+vz[i]*vz[i];
			}
			kin=0,5*m*(svxm+svym+svzm); //kinetic energy
			Tf=2*kin/(3*boltz*N); //final temperature
		}
	}
}

/*INPUTS:dim=numb. atoms; x,y,z=atom positions (arrays); c=eigenvectors (matrix with each column as the n-th eigenvector); rc=cut-off radius;
nnear=number of nearest neighbours (nn) to i-th atom (array); inear=label of j-th nn to i-th atom (matrix); fx,fy,fz forces on each atom (arrays);
maxnn=max number of nn */
int forces(int dim,int norbs,double *x,double *y,double *z,double *c,double rc,int *nnear,int *inear,double *fx,double *fy,double *fz)
{ int k,i,j,l,lp,n,m; /* dummy indeces for cycles*/
  double dd[3],ddm,ddmrx,ddmrrx,ddmlx,ddmllx,ddmry,ddmrry,ddmly,ddmlly,ddmrz,ddmrrz,ddmlz,ddmllz,ddrx[3],ddrrx[3],ddlx[3],ddllx[3],ddry[3],ddrry[3],ddly[3],ddlly[3],ddrz[3],ddrrz[3],ddlz[3],ddllz[3],drepx,drepy,drepz,h=rc/1000,sumphinn,sumphi;

  
  for(i=0;i<dim;i++){ /*initialisation of forces*/
    fx[i]=0;
    fy[i]=0;
    fz[i]=0;
    sumphinn=0;
  }

  for(i=0;i<dim;i++){ /*Cycle to compute band structure forces on atom i*/
    sumphi=0;
    	for(k=0;k<nnear[i];k++){
	dd[1]=abs(x[i]-x[inear[i*dim+k]]); /*Definition of vector distances*/
	dd[2]=abs(y[i]-y[inear[i*dim+k]]);
	dd[3]=abs(z[i]-z[inear[i*dim+k]]);

	ddm=sqrt(dd[1]*dd[1]+dd[2]*dd[2]+dd[2]*dd[2]); /*Modulus of distance*/
	
	sumphi=sumphi+o(ddm);
	}
    for(j=0;j<nnear[i];j++){ /*Cycle spanning the nearest neighbours of i*/
      if(j!=i){ /*Check to avoid self interaction (redundant, but saves operations)*/
	dd[1]=abs(x[i]-x[inear[i*dim+j]]); /*Definition of vector distances*/
	dd[2]=abs(y[i]-y[inear[i*dim+j]]);
	dd[3]=abs(z[i]-z[inear[i*dim+j]]);

	ddm=sqrt(dd[1]*dd[1]+dd[2]*dd[2]+dd[2]*dd[2]); /*Modulus of distance*/
	
	for(k=0;k<3;k++){ /*Initialisation of vector distances to perform derivatives */
	  ddrx[k]=dd[k]; 
	  ddrrx[k]=dd[k];
	  ddlx[k]=dd[k];
	  ddllx[k]=dd[k];

	  ddry[k]=dd[k];
	  ddrry[k]=dd[k];
	  ddly[k]=dd[k];
	  ddlly[k]=dd[k];
	
	  ddrz[k]=dd[k];
	  ddrrz[k]=dd[k];
	  ddlz[k]=dd[k];
	  ddllz[k]=dd[k];
	}
	
	ddrx[1]=abs(x[i]+h-x[j]);  /*Definition of variables to take derivatives (r->x+h,rr->x+2h,l->x-h,ll->x-h)*/
	ddrrx[1]=abs(x[i]+2*h-x[j]);
	ddlx[1]=abs(x[i]-h-x[j]);
	ddllx[1]=abs(x[i]-2*h-x[j]);

	ddry[2]=abs(y[i]+h-y[j]);
	ddrry[2]=abs(y[i]+2*h-y[j]);
	ddly[2]=abs(y[i]-h-y[j]);
	ddlly[2]=abs(y[i]-2*h-y[j]);

	ddrz[3]=abs(z[i]+h-z[j]);
	ddrrz[3]=abs(z[i]+2*h-z[j]);
	ddlz[3]=abs(z[i]-h-z[j]);
	ddllz[3]=abs(z[i]-2*h-z[j]);

	ddmrx=sqrt(ddrx[1]*ddrx[1]+ddrx[2]*ddrx[2]+ddrx[2]*ddrx[2]);
	ddmry=sqrt(ddry[1]*ddry[1]+ddry[2]*ddry[2]+ddry[2]*ddry[2]);
	ddmrz=sqrt(ddrx[1]*ddrz[1]+ddrz[2]*ddrz[2]+ddrz[2]*ddrz[2]);

	ddmrrx=sqrt(ddrrx[1]*ddrrx[1]+ddrrx[2]*ddrrx[2]+ddrrx[2]*ddrrx[2]);
	ddmrry=sqrt(ddrry[1]*ddrry[1]+ddrry[2]*ddrry[2]+ddrry[2]*ddrry[2]);
	ddmrrz=sqrt(ddrrz[1]*ddrrz[1]+ddrrz[2]*ddrrz[2]+ddrrz[2]*ddrrz[2]);
      
	ddmlx=sqrt(ddlx[1]*ddlx[1]+ddlx[2]*ddlx[2]+ddlx[2]*ddlx[2]);
	ddmly=sqrt(ddly[1]*ddly[1]+ddly[2]*ddly[2]+ddly[2]*ddly[2]);
	ddmlz=sqrt(ddlz[1]*ddlz[1]+ddlz[2]*ddlz[2]+ddlz[2]*ddlz[2]);
      
	ddmllx=sqrt(ddllx[1]*ddllx[1]+ddllx[2]*ddllx[2]+ddllx[2]*ddllx[2]);
	ddmlly=sqrt(ddlly[1]*ddlly[1]+ddlly[2]*ddlly[2]+ddlly[2]*ddlly[2]);
	ddmllz=sqrt(ddllz[1]*ddllz[1]+ddllz[2]*ddllz[2]+ddllz[2]*ddllz[2]);

	for(l=0;l<norbs;l++){ /*Cycle spanning the first orbital type*/
	  for(lp=0;l<norbs;l++){ /*Cycle spanning the second orbital type*/
	    for(n=0;n<dim;n++){ /*Cycle spanning the level of the level of the eigenvector*/
	      fx[i]=fx[i]-2*(-Gethijab(i,j,l,lp,ddrrx,6,6)+8*Gethijab(i,j,l,lp,ddrx,6,6)-8*Gethijab(i,j,l,lp,ddlx,6,6)+Gethijab(i,j,l,lp,ddllx,6,6))/(12*h)
		*c[l*dim+n]*c[lp*dim+n];
	      fy[i]=fy[i]-2*(-Gethijab(i,j,l,lp,ddrry,6,6)+8*Gethijab(i,j,l,lp,ddry,6,6)-8*Gethijab(i,j,l,lp,ddly,6,6)+Gethijab(i,j,l,lp,ddlly,6,6))/(12*h)
		*c[l*dim+n]*c[lp*dim+n];
	      fz[i]=fz[i]-2*(-Gethijab(i,j,l,lp,ddrrz,6,6)+8*Gethijab(i,j,l,lp,ddrz,6,6)-8*Gethijab(i,j,l,lp,ddlz,6,6)+Gethijab(i,j,l,lp,ddllz,6,6))/(12*h)
		*c[l*dim+n]*c[lp*dim+n];
	    }
	  }
	}

	sumphinn=0;

	for(m=0;m<nnear[inear[i*dim+j]];m++){
	    dd[1]=abs(x[inear[i*dim+j]]-x[inear[inear[i*dim+j]*dim+m]]); /*Definition of vector distances*/
	    dd[2]=abs(y[inear[i*dim+j]]-y[inear[inear[i*dim+j]*dim+m]]);
	    dd[3]=abs(z[inear[i*dim+j]]-z[inear[inear[i*dim+j]*dim+m]]);
		    
	    ddm=sqrt(dd[1]*dd[1]+dd[2]*dd[2]+dd[2]*dd[2]); /*Modulus of distance*/
	    
	    sumphinn=sumphinn+o(ddm);
	}
	      /*calculation of repuslve forces*/
	fx[i]=fx[i]-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrrx)+8*o(ddmrx)-8*o(ddmlx)+o(ddmllx))/(12*h);
	fy[i]=fy[i]-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrry)+8*o(ddmry)-8*o(ddmly)+o(ddmlly))/(12*h);
	fz[i]=fz[i]-(d_f0(sumphinn)+d_f0(sumphi))*(-o(ddmrrz)+8*o(ddmrz)-8*o(ddmlz)+o(ddmllz))/(12*h);
	
	
      }
    }
  }

  return 0;

}

void near_neigh(int N, double *x, double *y, double *z, double rc, int *nnear, int *inear, double sx, double sy, double sz)
{  //determine the nearest neighbours for each atom
	double dx,dy,dz,dist;
	for (int i=0; i<N; i++) { nnear[i]=0; }
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			dx=x[i]-x[j];
			dx=dx-sx*round(dx/sx);//if dx>sx, the distance to N.N. is dx-sx
			dy=y[i]-y[j];
			dy=dy-sy*round(dy/sy);
			dz=z[i]-z[j];
			dz=dz-sz*round(dz/sz);
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			if (dist<rc && i!=j) //add the atom j to the nearest neighbour list of i if this holds
			{
				nnear[i]=nnear[i]+1;
				inear[i*N+nnear[i]]=j; //a matrix with i rows, nnear[i] (no of nearest neighbours) columns
			}
		}
	}
}

void velocity(int N, double *mass, double *vx, double *vy, double *vz, double T)
{  //calculate velocities
	double c, ke, m=mass[1],boltz=1/11603;//1.38*pow(10,-23); //for now, each mass is the same
	double vxtot=0,vytot=0,vztot=0,msvx=0,msvy=0,msvz=0,vxavg,vyavg,vzavg,Tp;
	double g=sqrt(3*boltz*T/m);
	for (int i=0; i<N; i++)
	{
	  
	  srand (1);
		vx[i]=0; vy[i]=0; vz[i]=0;
//		c[i]=sqrt(3*boltz*T/m[i])
//		vx[i]=c[i]*(2*rand(0)-1); //random number generator required
//		vy[i]=c[i]*(2*rand(0)-1);
//		vz[i]=c[i]*(2*rand(0)-1);
		vx[i]=c*(2*rand()-1);
		vy[i]=c*(2*rand()-1);
		vz[i]=c*(2*rand()-1);
		vxtot=vxtot+vx[i];
		vytot=vytot+vy[i];
		vztot=vztot+vz[i];
	}
	vxavg=vxtot/N;
	vyavg=vytot/N;
	vzavg=vztot/N;

	for (int i=0; i<N; i++)
	{
		vx[i]=vx[i]-vxavg;
		vy[i]=vy[i]-vyavg;
		vz[i]=vz[i]-vzavg;
		msvx=msvx+vx[i]*vx[i];//mean square velocity
		msvy=msvy+vy[i]*vy[i];
		msvz=msvz+vz[i]*vz[i];
	}
	ke=0.5*m*(msvx+msvy+msvz); //must be adapted for different masses
	Tp=2*ke/(3*boltz);
	for(int i=0; i<N; i++)
	{
		vx[i]=vx[i]*sqrt(T/Tp);
		vy[i]=vy[i]*sqrt(T/Tp);
		vz[i]=vz[i]*sqrt(T/Tp);
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

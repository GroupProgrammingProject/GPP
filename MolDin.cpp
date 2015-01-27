#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <cstdlib> 
#include <stdio.h>
using namespace std;

int forces(int dim,int norbs,double *x,double *y,double *z,double H*,double c*,double rc,int *nnear,int *inear)
{ int i,j,l,lp,n;
  double dd,ddx,ddy,ddz,dsx,dsy,dsz,dphix,dphiy,dphiz,dphi;
  
  for(i=0;i<dim:i++){
    fx[i]=0;
    fy[i]=0;
    fz[i]=0;
  }
  for(i=0;i<dim;i++){
    for(j=0;j<nnear[i];j++){
      ddx=pow(x[i]-x[inear[i][j]],2);
      ddy=pow(y[i]-y[inear[i][j]],2);
      ddz=pow(z[i]-z[inear[i][j]],2);
      dd=sqrt(ddx+ddy+ddz);

      dsx=Derivsx(dd);
      dsy=Derivsy(dd);
      dsz=Derivsz(dd);

      dphix=Derivphix(dd);
      dphiy=Derivphiy(dd);
      dphiz=Derivphiz(dd);

      for(l=0;l<norbs;l++){
	for(lp=0;l<norbs;l++){
	  for(n=0;n<dim;n++){
	    fx[i]=fx[i]-2*H[i+l][j+lp]*dsx*c[l][n]*c[lp][n];
	    fy[i]=fy[i]-2*H[i+l][j+lp]*dsy*c[l][n]*c[lp][n];
	    fz[i]=fz[i]-2*H[i+l][j+lp]*dsz*c[l][n]*c[lp][n];
	  }
	}
      }
      fx[i]=fx[i]-dphix;
      fy[i]=fy[i]-dphiy;
      fz[i]=fz[i]-dphiz;
    }
  }

  return 0;

void near_neigh(int N, double *x, double *y, double *z, double rc, double *nnear, double *inear, double sx, double sy, double sz)
{
	double dist;
	for (int i=0; i<N; i++) { nnear[i]=0; }
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			dx=x[i]-x[j];
			dx=dx-sx*round(dx/sx);
			dy=y[i]-y[j];
			dy=dy-sy*round(dy/sy);
			dz=z[i]-z[j];
			dz=dz-sz*round(dz/sz);
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			if (dist<rc && i!=j)
			{
				nnear[i]=nnear[i]+1;
				inear[i][nnear[i]]=j;
			}
		}
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

#include "../include/chain.h"
//Nose-Hoover 2-thermostat chain 

void chain(int N,double m,double T,double dt,double q1,double q2,double &xi1,double &xi2,double &vxi1,double &vxi2,std::vector<double>* vx,std::vector<double>* vy,std::vector<double>* vz,double &kin)
{ //vxi1,vxi2,xi1,xi2 are coordinates introduced for the Nose-Hoover thermostat
	double boltz=1./11603,g1,g2;
	g2=(q1*vxi1*vxi1-boltz*T)/q2;
	vxi2=vxi2+g2*dt/4; //perform operations that change each term in succession (see Frenkel+Smit, Appendix E.2)
	vxi1=vxi1*exp(-vxi2*dt/8);
	g1=(2*kin-3*N*T)/q1;
	vxi1=vxi1+g1*dt/4;
	vxi1=vxi1*exp(-vxi2*dt/8);
	xi1=xi1-vxi1*dt/2;
	xi2=xi2-vxi2*dt/2;
	s=exp(-vxi1*dt/2);//scaling factor for velocities
	for(int i=0; i<N; i++)
	{
		(*vx).at(i)=s*(*vx).at(i);
		(*vy).at(i)=s*(*vy).at(i);
		(*vz).at(i)=s*(*vz).at(i);
	}
	kin=kin*s*s; //rescale kinetic energy
	vxi1=vxi1*exp(-vxi2*dt/8);
	g1=(2*kin-3*N*T)/q1;
	vxi1=vxi1+g1*dt/4;
	vxi1=vxi1*exp(-vxi2*dt/8);
	g2=(q1*vxi1*vxi1-boltz*T)/q2;
	vxi2=vxi2+g2*dt/4;
}

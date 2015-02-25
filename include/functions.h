#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <iostream>
#include <vector>
//scaling functions, notations as in the paper 

//tail function for s(r)
double ts (double r){
  double c0=6.7392620074314*pow(10,-3);
  double c1=-8.1885359517898*pow(10,-2);
  double c2=0.1932365259144;
  double c3=0.3542874332380;
  double r1=2.45;

  double ts=pow((r-r1),3)*c3+pow((r-r1),2)*c2+(r-r1)*c1+c0;
  return ts;
}

//derivative of s(r) with respect to dx
double dts (double r,double dx){
  double c0=6.7392620074314*pow(10,-3);
  double c1=-8.1885359517898*pow(10,-2);
  double c2=0.1932365259144;
  double c3=0.3542874332380;
  double r1=2.45;

  double dts=dx*(c1+2*c2*(r-r1)+3*c3*(r-r1)*(r-r1))/r;
  return dts;
}

//tail function for phi(r)
double to (double r){
  double c0=2.2504290109*pow(10,-8);
  double c1=-1.4408640561*pow(10,-6);
  double c2=2.1043303374*pow(10,-5);
  double c3=6.6024390226*pow(10,-5);
  double d1=2.57;

  double to=pow((r-d1),3)*c3+pow((r-d1),2)*c2+(r-d1)*c1+c0;
  return to;
}

//derivative of to
double dto (double r, double dx){
  double c0=2.2504290109*pow(10,-8);
  double c1=-1.4408640561*pow(10,-6);
  double c2=2.1043303374*pow(10,-5);
  double c3=6.6024390226*pow(10,-5);
  double d1=2.57;

  double dto=dx*(c1+2*c2*(r-d1)+3*c3*(r-d1)*(r-d1))/r;
  return dto;
}
 
//function s(r), including the tail and the atom types
double s (double r){
  double r0=1.536329; //nearest-neighbour atomic separation
  double n=2;
  double nc=6.5;
  double rc=2.18;
  double r1=2.45;
  double rcut=2.6;
  double S;
  //constants 
    if (r<r1){
      S=pow(r0/r,n)*exp(n*(-pow(r/rc,nc)+pow(r0/rc,nc)));
    }
    else if(r>=r1 && r<rcut){
      S=ts(r);
    }
   else{
    return 0;
   }
	return S;
}

//derivative of s(r) with respect to dx
double ds (double r,double dx){
  double r0=1.536329; //nearest-neighbour atomic separation
  double n=2;
  double nc=6.5;
  double rc=2.18;
  double r1=2.45;
  double rcut=2.6;
  double dS;
  //constants 
    if (r<r1){
      dS=-n*s(r)*(1+nc*pow(r/rc,nc))*dx/(r*r);
    }
    else if(r>=r1 && r<rcut){
      dS=dts(r,dx);
    }
   else{
    return 0;
   }
	return dS;
}

//phi(r) functions including the tail function
double o (double r){
  double phi0=8.18555;
  double m=3.30304;
  double mc=8.6655;
  double dc=2.1052;
  double d0=1.64;
  double d1=2.57;
  double dcut=2.6;
  double phi;

  if (r<d1){
    phi=phi0*pow(d0/r,m)*exp(m*(-pow(r/dc,mc)+pow(d0/dc,mc)));
  }
  else if(r>=d1 && r<dcut){
    phi=to(r);
  }
  else{
  		phi=0;
  }
  return phi;
}

//f polynomial needed for repulsive energy
double f0(double x){
  double c0=-2.5909765118191;
  double c1=0.5721151498619;
  double c2=-1.7896349903996*pow(10,-3);
  double c3=2.3539221516757*pow(10,-5);
  double c4=-1.24251169551587*pow(10,-7);

  double f0=pow(x,4)*c4+pow(x,3)*c3+pow(x,2)*c2+(x)*c1+c0;
  return f0;
}

//derevative of f polynomial needed for MD
double d_f0(double x){
  double c1=0.5721151498619;
  double c2=-1.7896349903996*pow(10,-3);
  double c3=2.3539221516757*pow(10,-5);
  double c4=-1.24251169551587*pow(10,-7);

  double f0=pow(x,3)*c4*4+pow(x,2)*c3*3+pow(x,1)*c2*2+c1;
  return f0;
}

//derivative of phi(r) needed for MD
double d_o(double r,double dx){
  double phi0=8.18555;
  double m=3.30304;
  double mc=8.6655;
  double dc=2.1052;
  double d0=1.64;
  double d1=2.57;
  double dcut=2.6;
  double dphi;

  if (r<d1){
    dphi=-(m/r)*o(r)*(1+mc*pow((r/dc),mc))*(dx/r);;
  }
  else if(r>=d1 && r<dcut){
    dphi=dto(r,dx);
  }
  else{
    dphi=0;
  }
  return dphi;
}

//X function - argument of f polynomial consisting SUM_over_j ( phi(r_ij))
double X (Eigen::MatrixXd* modr,int n, int i ){
  double r;
  double x=0;
  for (int j=0;j<n;j++){
    if(i!=j){
      r=(*modr)(i,j);
      if(r>1e-5){x=x+o(r);}
    };
  };
  return x;
};

//total repulsive energy - sum over all atoms of f(X_i)
double Erep (Eigen::MatrixXd* modr){
  int N=sqrt(modr->size());
  double x;
  double total=0.0;
  for(int i=0;i<N;i++){
    x=X(modr,N,i);
    total=total+f0(x);
  }
  return total;
}


#endif

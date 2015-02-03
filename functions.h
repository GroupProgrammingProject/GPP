#ifdef functions_h
#define functions_h

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
 
//function s(r), including the tail and the atom types
double s (double r,int type1, int type2){
  double r0=1.536329; //nearest-neighbour atomic separation
  double n=2;
  double nc=6.5;
  double rc=2.18;
  double r1=2.45;
  double S;
  //constants 
   if(type1==6 &&type2==6){
    if (r<r1){
      S=pow(r0/r,n)*exp(n*(-pow(r/rc,nc)+pow(r0/rc,nc)));
    }
    else{
      S=ts(r);
    }
    }
 
   else{
    return 0;
  }
  return S;
}

//phi(r) functions including the til function
double o (double r){
  double o0=8.18555;
  double m=3.30304;
  double mc=8.6655;
  double dc=2.1052;
  double d0=1.64;
  double d1=2.57;
  double O;

  if (r<d1){
    O=o0*pow(d0/r,m)*exp(m*(-pow(r/dc,mc)+pow(d0/dc,mc)));
  }
  else{
    O=to(r);
  }
  return O;
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

//X function - argument of f polynomial consisting SUM_over_j ( phi(r_ij))
double X (std::vector<int> &type, std::vector<double> &rx, std::vector<double> &ry,std::vector<double> &rz, int i ){
  double r;
  double x=0;
  int N=rx.size();
  std::cout<<" N "<<N<<std::endl;
  for (int j=0;j<N;j++){
    if(i!=j){
      r=sqrt(pow((rx[i]-rx[j]),2)+pow((ry[i]-ry[j]),2)+pow((rz[i]-rz[j]),2));
      std::cout<<" I: "<<i<<" J: "<<j<<" r: "<<r<<std::endl;
      x=x+o(r);
    };
  };
  return x;
};

//total repulsive energy - sum over all atoms of f(X_i)
double Erep (std::vector<int> &type, std::vector<double> &rx, std::vector<double> &ry,std::vector<double> &rz){
  int N=rx.size();
  std::cout<<" int N "<<N<<std::endl;
  double x;
  double total;
  for(int i=0;i<N;i++){
    std::cout<<"call X "<<std::endl;
    x=X(type,rx,ry,rz,i);
    total=total+f0(x);
  }
  return total;
}

//some test functions, diregard them 
void test (std::vector<double> &wektor){
  int N=wektor.size();
  for (int i=0;i<N;i++){
    std::cout<<wektor[i]<<std::endl;
  }

}


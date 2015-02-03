double ts (double r){
  double c0=6.7392620074314*pow(10,-3);
  double c1=-8.1885359517898*pow(10,-2);
  double c2=0.1932365259144;
  double c3=0.3542874332380;
  double r1=2.45;

  double ts=pow((r-r1),3)*c3+pow((r-r1),2)*c2+(r-r1)*c1+c0;
  return ts;
}

double S(double r, int typei, int typej) {
  double r0=1.536329; //nearest-neighbour atomic separation
  double n=2;
  double nc=6.5;
  double rc=2.18;
  double r1=2.45;
  double S;
  //constants 

  if (r<r1){
    S=pow(r0/r,n)*exp(n*(-pow(r/rc,nc)+pow(r0/rc,nc)));
  }
  else{
    S=ts(r);
  }
  return S;
}
